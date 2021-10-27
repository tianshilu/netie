netie=function(input_one_patient,sigma_square,alpha,beta,sigma_p_sqr,sigma_a_sqr=NULL,max_iter,multi_sample=FALSE,
               cellular_clock='variant_allele_frequency',
               cellular_prevalence_min=0.02,
               keep_mutations_number=2,
               keep_neoantigen_encoding_mutations_number=1){
  #cellular prevalence filter
  input_one_patient=input_one_patient[input_one_patient$cellular_prevalence>=cellular_prevalence_range[1] ,]
  if(all(input_one_patient$neo_load[!is.na(input_one_patient$cluster_id)]==0)){
    return(NA)
  }
  input_one_patient=input_one_patient[!is.na(input_one_patient$cluster_id),]
  
  #change cellular prevalence or vaf to cellular clock
  if(cellular_clock=='variant_allele_frequency'){
    colnames(input_one_patient)[colnames(input_one_patient)=='variant_allele_frequency']='cellular_clock'
  }else{
    colnames(input_one_patient)[colnames(input_one_patient)=='cellular_prevalence']='cellular_clock'
  }
  
  #multi_sample
  if(multi_sample==T){
    #same mutations have same neoantigens
    mutations=unlist(sapply(input_one_patient$mutation_id,function(x) paste(strsplit(x,' ')[[1]][2],
                                                                            strsplit(x,' ')[[1]][3])))
    input_one_patient$neo_load=unlist(sapply(mutations,function(x) max(input_one_patient[mutations==x,'neo_load'])))
    #find similar clones
    phi='1'
    clones=list()
    clones[[id='1']]=mutations[paste(input_one_patient$sample_id,input_one_patient$cluster_id)==
                                 paste(input_one_patient$sample_id,input_one_patient$cluster_id)[1]]
    
    for(each_clone in unique(paste(input_one_patient$sample_id,input_one_patient$cluster_id))[-1]){
      mutations_one_clone=mutations[paste(input_one_patient$sample_id,input_one_patient$cluster_id)==each_clone]
      phi_tmp=unlist(sapply(1:length(clones),function(x) {uniq_clone=clones[[x]]
      shared_mutations=intersect(uniq_clone,mutations_one_clone)
      #if shared mutations are 50% or more than considered as same clone
      if(length(shared_mutations)/length(uniq_clone)>0.5 & 
         length(shared_mutations)/length(mutations_one_clone)>0.5){
        return(names(clones)[x])
      }
      }),use.names = F)
      if(!is.null(phi_tmp)){
        phi=c(phi,phi_tmp)
      }else{
        phi_tmp=max(as.numeric(names(clones)))+1
        phi=c(phi,phi_tmp)
        clones[[id=as.character(phi_tmp)]]=mutations_one_clone
      }
    }
    names(phi)=unique(paste(input_one_patient$sample_id,input_one_patient$cluster_id))
  }
  
  if(length(unique(input_one_patient$cluster_id))>1){
    if(is.null(sigma_a_sqr)){
      non_zero_neo_avg=sapply(unique(input_one_patient$cluster_id),function(x) 
        mean(input_one_patient[input_one_patient$cluster_id==x 
                               & input_one_patient$neo_load!=0,'neo_load']))
      non_zero_neo_avg[is.nan(non_zero_neo_avg)]=0
      sigma_a_sqr=sd(log(non_zero_neo_avg+1))^2*10
      if(sigma_a_sqr==0){
        sigma_a_sqr=1
      }
    }}else{
      sigma_a_sqr=1 
    }
  #check sigma_square >> sigma_a_sqr
  if(sigma_square<100*sigma_a_sqr){
    print("sigma square should be much more larger than sigma a square!")
    stop()
  }
  #alpha should be larger than beta
  if(alpha<=beta){
    print("alpha should be larger than beta!")
    stop()
  }
  
  if(multi_sample==T){
    #keep mutations with vaf>0.5 in any samples
    max_vaf=unlist(sapply(mutations,function(x) max(input_one_patient[mutations==x,'variant_allele_frequency'])))
    input_one_patient=input_one_patient[max_vaf>0.05,]
  }else{
    #only keep mutations with vaf>0.05 in single samples
    input_one_patient=input_one_patient[input_one_patient$variant_allele_frequency>0.05,]
  }
  #keep clusters with a customized number of mutations with neoantigens
  tmp=table(input_one_patient$cluster_id[input_one_patient$neo_load>0]) 
  tmp=names(tmp[tmp>= keep_neoantigen_encoding_mutations_number])
  input_one_patient=input_one_patient[input_one_patient$cluster_id %in% tmp,]
  #only keep clones with >=2 mutations
  tmp=table(input_one_patient$cluster_id)
  tmp=names(tmp[tmp>=keep_mutations_number])
  input_one_patient=input_one_patient[input_one_patient$cluster_id %in% tmp,]
  if (dim(input_one_patient)[1]==0) {return(NA)} 
  
  if(multi_sample==T){
    input_one_patient$phi=as.numeric(phi[paste(input_one_patient$sample_id,input_one_patient$cluster_id)])
    input_one_patient$cluster_id=as.numeric(factor(paste(input_one_patient$sample_id,
                                                         input_one_patient$cluster_id)))
    #map cluster with phi
    phi_cluster=input_one_patient[,c('cluster_id','phi')]
    phi_cluster=phi_cluster[!duplicated(phi_cluster$cluster_id),]
    rownames(phi_cluster)=as.character(phi_cluster$cluster_id)
    phi_cluster=phi_cluster[as.character(unique(input_one_patient$cluster_id)),]
  }else{
    input_one_patient$cluster_id=
      as.numeric(factor(input_one_patient$cluster_id))
  }
  input_one_patient[input_one_patient$neo_load>150,'neo_load']=150
  ########  initializaton  #############
  ac=bc=rep(0,length(unique(input_one_patient$cluster_id))) 
  pi=0.5
  a=0
  ########  iterations  ############
  zck_list=list()
  ac_list=list()
  bc_list=list()
  
  acp_rate_ac_list=list()
  acp_rate_bc_list=list()
  
  a_all=c()
  pi_all=c()
  
  for (iter in 1:max_iter)
  {
    if(iter/1000==round(iter/1000)){
      cat(paste("Iteration",iter,"\n"))
    }
    #record acceptance rate
    acp_rate_ac=rep(FALSE,length(unique(input_one_patient$cluster_id)))
    acp_rate_bc=rep(FALSE,length(unique(input_one_patient$cluster_id)))
    
    zck_df=input_one_patient[,c('mutation_id','cluster_id')];zck_df$zck=1
    if(multi_sample==T){
      for(p in 1:length(unique(input_one_patient$phi))){
        input_each_phi=input_one_patient[input_one_patient$phi==unique(input_one_patient$phi)[p],]
        for(c in unique(input_each_phi$cluster_id)){
          #c:cluster_id
          input_each_clone=input_each_phi[input_each_phi$cluster_id==c,]
          
          vck=input_each_clone$variant_allele_frequency
          lambda=exp(ac[c]*vck+bc[c])
          nck=input_each_clone$neo_load
          
          #update zck
          r_tmp=pi*(nck==0)/(pi*(nck==0)+(1-pi)*dpois(nck,lambda,log=F))
          r_tmp_deno=pi*(nck==0)+(1-pi)*dpois(nck,lambda,log=F)
          r_tmp[r_tmp_deno==0]=0
          zck=1*(runif(length(nck),0,1)>r_tmp)
          names(zck)=input_each_clone$mutation_id
          zck_df$zck[zck_df$mutation_id %in% names(zck)]=zck
          
          #update bc
          bc_prim=rnorm(1,bc[c],sqrt(sigma_p_sqr))
          lambda_prim_b=exp(ac[c]*vck+bc_prim)
          lambda=exp(ac[c]*vck+bc[c])
          
          tmp_prim=sum((zck==1)*dpois(nck,lambda_prim_b,log = T))
          tmp=sum((zck==1)*dpois(nck,lambda,log = T))
          llhr_b=exp(tmp_prim-bc_prim^2/(2*sigma_square)-tmp+bc[c]^2/(2*sigma_square))
          
          acceptance_function_b=min(1,llhr_b) 
          
          u=runif(1,0,1)
          if(u<=acceptance_function_b){
            bc[c]=bc_prim
            acp_rate_bc[c]=TRUE
          }
        }
        input_each_phi$bc=bc[input_each_phi$cluster_id]
        input_each_phi$ac=ac[c]
        
        vck_phi=input_each_phi$variant_allele_frequency
        lambda_phi=exp(input_each_phi$ac*vck_phi+input_each_phi$bc)
        nck_phi=input_each_phi$neo_load
        
        zck_phi=zck_df[input_each_phi$mutation_id,'zck']
        
        #update ac
        ac_prim=rnorm(1,ac[c],sqrt(sigma_p_sqr))
        lambda_prim_a=exp(ac_prim*vck_phi+input_each_phi$bc)
        
        #calculate likelihood ratio for new ac and old ac  
        tmp_prim=sum((zck_phi==1)*dpois(nck_phi,lambda_prim_a,log = T))
        tmp=sum((zck_phi==1)*dpois(nck_phi,lambda_phi,log = T))
        
        
        if(length(table(input_one_patient$cluster_id))==1){
          #the patient only has one clone
          llhr_a=exp(tmp_prim-ac_prim^2/(2*sigma_square)-tmp+ac[c]^2/(2*sigma_square))
        }else{
          llhr_a=exp(tmp_prim-(ac_prim-a)^2/(2*sigma_a_sqr)-tmp+(ac[c]-a)^2/(2*sigma_a_sqr))
        }
        
        acceptance_function_a=min(1,llhr_a)
        
        u=runif(1,0,1)
        if(u<=acceptance_function_a){
          ac[phi_cluster$phi==unique(input_each_clone$phi)]=ac_prim
          acp_rate_ac[c]=TRUE
        }
      }
      #update pi
      pi=rbeta(1,alpha+sum((zck_df$zck==0)*(input_one_patient$neo_load==0)),beta+sum(zck_df$zck==1))
      
      #update a
      A=1/sigma_square+length(unique(input_one_patient$phi))/sigma_a_sqr
      B=sum(ac[!duplicated(phi_cluster$phi)])/sigma_a_sqr
      
      a=rnorm(1,B/A,sqrt(1/A))
      
      #save results
      ac_list[[iter]]=ac
      bc_list[[iter]]=bc
      zck_list[[iter]]=zck_df$zck
      acp_rate_ac_list[[iter]]=acp_rate_ac
      acp_rate_bc_list[[iter]]=acp_rate_bc
      
      a_all=c(a_all,a)
      pi_all=c(pi_all,pi)
    }else{
      for(c in 1:length(unique(input_one_patient$cluster_id))){
        input_each_clone=input_one_patient[input_one_patient$cluster_id==unique(input_one_patient$cluster_id)[c],]
        
        vck=input_each_clone$variant_allele_frequency
        
        lambda=exp(ac[c]*vck+bc[c])
        nck=input_each_clone$neo_load
        
        #update zck
        r_tmp=pi*(nck==0)/(pi*(nck==0)+(1-pi)*dpois(nck,lambda,log=F))
        r_tmp_deno=pi*(nck==0)+(1-pi)*dpois(nck,lambda,log=F)
        r_tmp[r_tmp_deno==0]=0
        zck=1*(runif(length(nck),0,1)>r_tmp)
        names(zck)=input_each_clone$mutation_id
        zck_df$zck[zck_df$mutation_id %in% names(zck)]=zck
        
        #update ac
        ac_prim=rnorm(1,ac[c],sqrt(sigma_p_sqr))
        lambda_prim_a=exp(ac_prim*vck+bc[c])
        
        #calculate likelihood ratio for new ac and old ac  
        tmp_prim=sum((zck==1)*dpois(nck,lambda_prim_a,log = T))
        tmp=sum((zck==1)*dpois(nck,lambda,log = T))
        
        if(length(table(input_one_patient$cluster_id))==1){
          #the patient only has one clone
          llhr_a=exp(tmp_prim-ac_prim^2/(2*sigma_square)-tmp+ac[c]^2/(2*sigma_square))
        }else{
          llhr_a=exp(tmp_prim-(ac_prim-a)^2/(2*sigma_a_sqr)-tmp+(ac[c]-a)^2/(2*sigma_a_sqr))
        }
        
        acceptance_function_a=min(1,llhr_a)
        
        u=runif(1,0,1)
        if(u<=acceptance_function_a){
          ac[c]=ac_prim
          acp_rate_ac[c]=TRUE
        }
        
        
        #update bc
        bc_prim=rnorm(1,bc[c],sqrt(sigma_p_sqr))
        lambda_prim_b=exp(ac[c]*vck+bc_prim)
        lambda=exp(ac[c]*vck+bc[c])
        
        tmp_prim=sum((zck==1)*dpois(nck,lambda_prim_b,log = T))
        tmp=sum((zck==1)*dpois(nck,lambda,log = T))
        llhr_b=exp(tmp_prim-bc_prim^2/(2*sigma_square)-tmp+bc[c]^2/(2*sigma_square))
        
        acceptance_function_b=min(1,llhr_b) 
        
        u=runif(1,0,1)
        if(u<=acceptance_function_b){
          bc[c]=bc_prim
          acp_rate_bc[c]=TRUE
        }
      }
      #update pi
      pi=rbeta(1,alpha+sum((zck_df$zck==0)*(input_one_patient$neo_load==0)),beta+sum(zck_df$zck==1))
      
      #update a
      A=1/sigma_square+length(unique(input_one_patient$cluster_id))/sigma_a_sqr
      B=sum(ac)/sigma_a_sqr
      
      a=rnorm(1,B/A,sqrt(1/A))
      
      #save results
      ac_list[[iter]]=ac
      bc_list[[iter]]=bc
      zck_list[[iter]]=zck_df$zck
      acp_rate_ac_list[[iter]]=acp_rate_ac
      acp_rate_bc_list[[iter]]=acp_rate_bc
      
      a_all=c(a_all,a)
      pi_all=c(pi_all,pi)
    }
  }
  #take average
  keep=round(max_iter/2):max_iter
  ac_final=Reduce("+",ac_list[keep])/length(keep)
  bc_final=Reduce("+",bc_list[keep])/length(keep)
  zck_df_final=round(Reduce("+",zck_list[keep])/length(keep))
  names(zck_df_final)=zck_df$mutation_id
  
  ac_rate=Reduce("+",acp_rate_ac_list[keep])/length(keep)
  bc_rate=Reduce("+",acp_rate_bc_list[keep])/length(keep)
  
  a_final=mean(a_all[keep])
  pi_final=mean(pi_all[keep])
  
  if(multi_sample==TRUE){
    final_parameters=list(ac=cbind(phi_cluster,ac_final),a=a_final)
  }else{
    final_parameters=list(ac=ac_final,a=a_final)
  }
  result=list('final_parameters'=final_parameters)
  return(result)
}
