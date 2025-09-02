## Introduction

The code in this repository can be used to reproduce the analysis in Morrison et al (2024), "Genetic Factor Analysis for Characterizing Phenome-Wide Patterns of Genetic Pleiotropy". To run your own GFA analysis, please go to [the main GFA repository[(https://github.com/jean997/GFA). You will find detailed instructions as well as a link to our Snakemake pipeline which will run all steps of the analysis including data harmonization. 

## Pre-Requisite Software
We used the following software to run the analysis in this directory. We have provided version numbers, however, it is most likely unnecessary to match our version numbers in order to reproduce our analysis. 

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- python 3.13.5
- R 4.4.3

R-packages: 
- VariantAnnotation (Bioconductor) 1.46.0
- gwasvcf (github, mrcieu) 0.1.4
- dplyr  (CRAN) 1.1.4
- readr (CRAN) 2.1.5
- purrr (CRAN) 1.1.0 
- stringr (CRAN) 1.5.1
- ieugwasr (github, mrcieu) 1.1.0
- MVMR (github WSpiller) 0.4.1
- viridis (CRAN)
- qqman (CRAN)
- topr (CRAN)
- GFA (github, jean997)  0.1.2.0333


## Reproduce the Applied analysis

### 1. Downloa data: 
Change directories to the data/ directory. 
Execute
```
./download_data.sh
```
This will download all of the necessary GWAS summary statistics from zenodo. It may take some time. 

### 2. Change directories to code/applied_analysis. 
We strongly recommend executing the snakemake pipeline on a compute cluster. 
This can be done by executing run-snakemake-gfa_v9.sh. However, you should make sure that 
the commands will work for your cluster setup first. Open this file in a text editor (e.g. Vim). 

If your cluster uses SLURM, you may not need to make any changes to the execution script. 
If you need to specify an account, you can do this on the --default-resources line, e.g.

```
--default-resources mem_mb=5000 runtime=120 account=myaccount \
```

If your cluster does not use slurm, you will need to use a different executor. See 
the Snakemake documentation for more details. 

### 3. Execute the snakemake pipeline using 

```
nohup ./run-snakemake-gfa_v9.sh &
```
Using nohup and the "&" will allow it to run uninterrupted as the process may take a few hours. 

4. Look at reults. Check the end of nohup.out to make sure the pipeline ran without errors. At the
end of the pipeline, you should have directories

```
data:
    gfa_intermediate_data
results:
    applied_analysis:
        gfa_results
        mr_results
        plots
```

The contents of the `results/applied_analysis/plots` directory should match the contents of the `expected_applied_result_plots` directory of this repository which has
been included for convenient checking. 
Fitted GFA objects are in the `gfa_results` subdirectory. The two sets of results presented in the manuscript are in 
`results/applied_analysis/gfa_results/metab_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS` and `results/applied_analysis/gfa_results/bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS`. All of the other results files were used for various supplementary comparisons. 
