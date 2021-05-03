# Netie
The Bayesian Hierarchical Model named Neoantgien-T cell interaction estimation (Netie) is to investigate the neoantigens observed in the patient tumors to esti- mate the history of the immune pressure on the evolution of the tumor clones. The estimation results will answer whether the host immune system has been conferring strong or weaker selection pressure on the tumor clones over the time of tumor development. This may give us a peak into the future of how the mutations and clones will evolve for that patient. The model is based on pyclone estimation results. Essentially, each clone is modelled separately, but sharing some random variables.
![preview](https://github.com/tianshilu/Netie/blob/main/flowchart.png)
## Dependencies
PyClone; R; 
## Guided Tutorial
Command: 
```
netie(input_data,sigma_square = 100000 ,
                                   alpha = 10,beta = 2,sigma_p_sqr = 0.1,sigma_a_sqr = NULL,max_iter =100000,multi_sample = T)
```
* input_data: a list with each data frame as the data for each patient. Each data frame consists 7 columns and each row is for one mutation. The 7 columns are mutation ID, sample ID, cluster ID, cellular prevalence, variant allele prevalence, variant allele frequency, and neoantigen load with column names as "mutation_id","sample_id","cluster_id","cellular_prevalence","variant_allele_frequency", and "neoantigen_load". 
* sigma_square                                  
