Netie
=====

An R package to infer the anti-tumor selection pressure for tumors.

Introduction
------------

The Bayesian Hierarchical Model named Neoantgien-T cell interaction
estimation (Netie) is to investigate the neoantigens observed in the
patient tumors to esti- mate the history of the immune pressure on the
evolution of the tumor clones. The estimation results will answer
whether the host immune system has been conferring strong or weaker
selection pressure on the tumor clones over the time of tumor
development. This may give us a peak into the future of how the
mutations and clones will evolve for that patient. The model is based on
pyclone estimation results. Essentially, each clone is modelled
separately, but sharing some random variables. 

![NetieFlowchart](https://github.com/tianshilu/Netie/blob/main/flowchart.jpg)

Please refer to our lab’s website for more information
<a href="https://qbrc.swmed.edu/labs/wanglab/software.php" class="uri">https://qbrc.swmed.edu/labs/wanglab/software.php</a>.

Installation of the package:
----------------------------

To install our package, you may simply execute the following codes.

``` r
# install.packages("devtools") 

devtools::install_github("tianshilu/netie", subdir = "netie") # don't forget to specify subdir!
```

Dependencies
------------

PyClone; R(version&gt;3.4)

Guided Tutorial
---------------

Command:

``` r
netie(input_data,sigma_square = 100000 ,
      alpha = 10,beta = 2,sigma_p_sqr = 0.1,sigma_a_sqr = NULL,max_iter =100000,
      cellular_clock='variant_allele_frequency',
      cellular_prevalence_min=0.02,
      keep_mutations_number=2,
      keep_neoantigen_encoding_mutations_number=1,
      multi_sample = T)
```

-   input\_data: a list with each data frame as the data for each
    patient. 
    
    Each data frame consists 7 columns and each row is for one
    mutation. The 7 columns are mutation ID, sample ID, cluster ID,
    cellular prevalence, variant allele prevalence, variant allele
    frequency, and neoantigen load with column names as
    “mutation\_id”,“sample\_id”,“cluster\_id”,“cellular\_prevalence”,“variant\_allele\_frequency”,
    and “neoantigen\_load”. 
    
    Please use PyClone or other softwares
    (<a href="https://github.com/tianshilu/Phylogenetic-Tree" class="uri">https://github.com/tianshilu/Phylogenetic-Tree</a>)
    to get information of cluster id and cellular prevalence (we recommend to keep mutations with sequencing depth more than 50 for clustering and netie inference). 
    
    Please use
    QBRC mutation calling pipeline
    (<a href="https://github.com/tianshilu/QBRC-Somatic-Pipeline" class="uri">https://github.com/tianshilu/QBRC-Somatic-Pipeline</a>)
    to call mutations for whole exome sequenicng; QBRC neoantigen
    calling pipeline
    (<a href="https://github.com/tianshilu/QBRC-Neoantigen-Pipeline" class="uri">https://github.com/tianshilu/QBRC-Neoantigen-Pipeline</a>)
    to call neoantigens for whole exome sequencing and RNA sequencing.

examples of input\_data:

``` r
example_input1= read.table('example_input1.txt', header=T,sep='\t',stringsAsFactors = F)
head(example_input1)
```

    ##                   mutation_id       sample_id cluster_id cellular_prevalence
    ## 1 TCGA-D3-A2JG-06 7 142881273 TCGA-D3-A2JG-06          0           0.3670081
    ## 2  TCGA-D3-A2JG-06 5 71016376 TCGA-D3-A2JG-06          0           0.3893199
    ## 3 TCGA-D3-A2JG-06 11 74015355 TCGA-D3-A2JG-06          0           0.3963911
    ## 4 TCGA-D3-A2JG-06 12 95914926 TCGA-D3-A2JG-06          0           0.3966120
    ## 5 TCGA-D3-A2JG-06 7 140453136 TCGA-D3-A2JG-06          0           0.3990956
    ## 6  TCGA-D3-A2JG-06 1 47746773 TCGA-D3-A2JG-06          0           0.3993479
    ##   cellular_prevalence_std variant_allele_frequency neo_load
    ## 1              0.07829626                0.1771429        0
    ## 2              0.09846122                0.2000000        0
    ## 3              0.11910942                0.2864583        0
    ## 4              0.12477982                0.3121387        5
    ## 5              0.12635598                0.3260870        1
    ## 6              0.11804737                0.2951807        0

``` r
library(ggplot2)

ggplot(example_input1,aes(variant_allele_frequency,neo_load))+ geom_point(colour = "seagreen3", size = 3) +labs(x = "Variant Allele Frequency",y="#Neoantigen per mutation")
```

![](https://github.com/tianshilu/Netie/blob/main/README_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
example_input2= read.table('example_input2.txt', header=T,sep='\t',stringsAsFactors = F)
head(example_input2)
```

    ##                         mutation_id       sample_id cluster_id
    ## 1 TCGA-BF-A3DN-01 GL000205.1 117561 TCGA-BF-A3DN-01          0
    ## 2       TCGA-BF-A3DN-01 1 153270485 TCGA-BF-A3DN-01          0
    ## 3         TCGA-BF-A3DN-01 2 1926514 TCGA-BF-A3DN-01          0
    ## 4        TCGA-BF-A3DN-01 4 47663875 TCGA-BF-A3DN-01          0
    ## 5       TCGA-BF-A3DN-01 17 11840675 TCGA-BF-A3DN-01          0
    ## 6        TCGA-BF-A3DN-01 7 93108821 TCGA-BF-A3DN-01          0
    ##   cellular_prevalence cellular_prevalence_std variant_allele_frequency neo_load
    ## 1           0.3441747              0.06237205                0.1290323        0
    ## 2           0.3450515              0.03324016                0.1538462        0
    ## 3           0.3470656              0.05007811                0.1372549        0
    ## 4           0.3529826              0.05362554                0.1600000        0
    ## 5           0.3537355              0.03360195                0.1655172        3
    ## 6           0.3588706              0.06354249                0.3085714        0

``` r
ggplot(example_input1,aes(variant_allele_frequency,neo_load))+ geom_point(colour = "blue", size = 3) +labs(x = "Variant Allele Frequency",y="#Neoantigen per mutation")
```

![](https://github.com/tianshilu/Netie/blob/main/README_files/figure-markdown_github/unnamed-chunk-3-2.png)

-   sigma\_square, alpha, beta, sigma\_p\_sqr, sigma\_a\_sqr:
    hyperparameters for prior distributions. Please refer to the paper
    for more details.

-   max\_iter: the iterations of Markov chain Monte Carlo.

-   cellular\_clock: choose to use cellular prevalence or variant allele frequency as the indicator of developmental time; variant allele frequency is the default.

-   cellular\_prevalence\_min: the minimal cutoff for cellular prevalence of clones; the default is 0.02.

-   keep\_mutations\_number: number of mutations you want to keep for inference; the minimal number is 2 and the default is 2.

-   keep\_neoantigen\_encoding\_mutations\_number: number of neoantigen-encoding mutations you want to keep for inference; the minimal number is 1 and the default is 1.

-   multi\_sample: use True if one patient has more than one sample.

Two example input datasets can be found here:
<a href="https://github.com/tianshilu/Netie/tree/main/example" class="uri">https://github.com/tianshilu/Netie/tree/main/example</a>

Output
------

The output is a list with the information of the anti-tumor selection
pressure for each clone ac and for the whole tumor a.

Two example output results can be found here:
<a href="https://github.com/tianshilu/Netie/tree/main/example" class="uri">https://github.com/tianshilu/Netie/tree/main/example</a>

examples of output

``` r
load('example_output1.RData')
print(example_output$a)
```

    ## [1] 1.577462

``` r
load('example_output2.RData')
print(example_output$a)
```

    ## [1] 3.885662
