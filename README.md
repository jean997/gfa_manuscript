## Introduction

The code in this repository can be used to reproduce the analysis in Morrison et al (2024), "Genetic Factor Analysis for Characterizing Phenome-Wide Patterns of Genetic Pleiotropy". To run your own GFA analysis, please go to [the main GFA repository[(https://github.com/jean997/GFA). You will find detailed instructions as well as a link to our Snakemake pipeline which will run all steps of the analysis including data harmonization. 

## Pre-Requisite Software
We used the following software to run the analysis in this directory. We have provided version numbers, however, it is most likely unnecessary to match our version numbers in order to reproduce our analysis. 

- [Snakemake](https://snakemake.readthedocs.io/en/stable/) 7.32.3 - used for applied anlaysis only
- [DSC](https://stephenslab.github.io/dsc-wiki/overview.html) 0.4.3.5 - used for simulations only
- python 3.11.5
- R 4.3.1

R-packages: 
- VariantAnnotation (Bioconductor) 1.46.0
- gwasvcf (github, mrcieu) 0.1.2
- dplyr  (CRAN) 1.1.4
- readr (CRAN) 2.1.5
- purrr (CRAN) 1.0.2 
- stringr (CRAN) 1.5.1
- ieugwasr (github, mrcieu) 0.1.5
- GFA (github, jean997)  0.1.2.0333

A full `sessionInfo()` output is below:
```
R version 4.3.1 (2023-06-16)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS:   /sw/pkgs/arc/stacks/gcc/10.3.0/R/4.3.1/lib64/R/lib/libRblas.so 
LAPACK: /sw/pkgs/arc/stacks/gcc/10.3.0/R/4.3.1/lib64/R/lib/libRlapack.so;  LAPACK version 3.11.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Detroit
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] gwasvcf_0.1.2               GFA_0.1.2.0333             
 [3] ieugwasr_0.1.5              stringr_1.5.1              
 [5] purrr_1.0.2                 readr_2.1.5                
 [7] rlang_1.1.2                 dplyr_1.1.4                
 [9] VariantAnnotation_1.46.0    Rsamtools_2.16.0           
[11] Biostrings_2.68.1           XVector_0.40.0             
[13] SummarizedExperiment_1.30.2 Biobase_2.60.0             
[15] GenomicRanges_1.52.0        GenomeInfoDb_1.36.2        
[17] IRanges_2.34.1              S4Vectors_0.38.1           
[19] MatrixGenerics_1.12.3       matrixStats_1.3.0          
[21] BiocGenerics_0.46.0        

loaded via a namespace (and not attached):
 [1] DBI_1.1.3                bitops_1.0-7             biomaRt_2.56.1          
 [4] magrittr_2.0.3           horseshoe_0.2.0          compiler_4.3.1          
 [7] RSQLite_2.3.1            GenomicFeatures_1.52.2   png_0.1-8               
[10] vctrs_0.6.4              reshape2_1.4.4           pkgconfig_2.0.3         
[13] crayon_1.5.3             fastmap_1.1.1            dbplyr_2.3.3            
[16] utf8_1.2.4               tzdb_0.4.0               bit_4.0.5               
[19] zlibbioc_1.46.0          cachem_1.0.8             trust_0.1-8             
[22] progress_1.2.3           blob_1.2.4               DelayedArray_0.26.7     
[25] BiocParallel_1.34.2      irlba_2.3.5.1            parallel_4.3.1          
[28] prettyunits_1.2.0        R6_2.5.1                 stringi_1.8.2           
[31] SQUAREM_2021.1           rtracklayer_1.60.1       Rcpp_1.0.12             
[34] Matrix_1.6-1             splines_4.3.1            tidyselect_1.2.1        
[37] abind_1.4-5              yaml_2.3.7               codetools_0.2-19        
[40] curl_5.0.2               lattice_0.21-8           tibble_3.2.1            
[43] plyr_1.8.9               KEGGREST_1.40.0          flashier_1.0.14         
[46] BiocFileCache_2.8.0      xml2_1.3.5               lpSolve_5.6.20          
[49] pillar_1.9.0             filelock_1.0.2           softImpute_1.4-1        
[52] generics_0.1.3           invgamma_1.1             RCurl_1.98-1.12         
[55] truncnorm_1.0-9          hms_1.1.3                ggplot2_3.5.1           
[58] munsell_0.5.1            scales_1.3.0             ashr_2.2-63             
[61] glue_1.7.0               tools_4.3.1              BiocIO_1.10.0           
[64] BSgenome_1.68.0          GenomicAlignments_1.36.0 XML_3.99-0.14           
[67] grid_4.3.1               tidyr_1.3.1              AnnotationDbi_1.62.2    
[70] colorspace_2.1-1         GenomeInfoDbData_1.2.10  deconvolveR_1.2-1       
[73] restfulr_0.0.15          cli_3.6.2                rappdirs_0.3.3          
[76] ebnm_1.1-34              fansi_1.0.6              mixsqp_0.3-54           
[79] S4Arrays_1.0.6           gtable_0.3.5             digest_0.6.33           
[82] rjson_0.2.21             memoise_2.0.1            lifecycle_1.0.4         
[85] httr_1.4.7               bit64_4.0.5             
```

## Requirements

The data used in the applied analysis occupies about 40G. The simulation experiments will produce almost 2 Tb of data. However, we have divided the simulations into chunks and most data is deleted after it is no longer needed. To fully reproduce the simulations, you will need to have at least X Gb of storage available. 

This analysis will run most quickly on a compute cluster. However, we have also provided instructions in the case that you do not have a compute cluster. Parallelization is most critical for the simulation portion of the analysis. For the applied analysis, parallelization speeds up data processing. Processing data for larger chromosomes takes about 30 minutes for the blood cell analysis and about 15 minutes for the metabolic traits analysis. If this data processing is done in sequence rather than in parallel, please allow about 10 hours for the analysis. If done in parallel, the whole analysis can be completed in under an hour.


## Instructions

First clone this repository to your computer. 

### Applied Analysis

1. Download the data. Change directories into the `data` directory. Use the following command to download the data from zenodo, unzip it and place it in the correct folders. Please note that this will download about 40Gb of data.
```
./download_data.sh
```

Most of the data size is the bloodcell trait data. If you would like to only replicate the metabolic traits analysis, instead run
```
./download_metab_only.sh
```



2. Change to the `code/applied_analysis` directory.
3. Use Snakemake to execute the analysis:
  3a. If you are using a compute cluster: Open the file `run-snakemake-gfa.sh` in a text editor (e.g. Vim). Make sure that the command given to the `--cluster` argument is appropriate for your cluster. If your cluster uses SLURM, you may not need to make any changes. However, if your cluster uses an account system, you may need to add `--account=[myaccount] \` to the cluster command. Run the following commands:
   3b. If you are not using a compute cluster: Run the following commands:

### Simulations

