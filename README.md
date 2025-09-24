## Introduction

The code in this repository can be used to reproduce the analysis in Morrison et al (2025), "Genetic Factor Analysis for Characterizing Phenome-Wide Patterns of Genetic Pleiotropy". To run your own GFA analysis, please go to [the main GFA repository](https://github.com/jean997/GFA). You will find detailed instructions as well as a link to our Snakemake pipeline which will run all steps of the analysis including data harmonization. 

## Pre-Requisite Software
You will need a conda installation and R. We used R 4.4.3 for this analysis. We have tested this repository in Linux only. If you want to run the whole 
analysis you will also need a compute cluster. The applied analysis requires downloading 42 Gb of data and will produce an additional 12 G of output files. 
The simulations will produce about 3.5T of data total. Most of this will be deleted automatically at the end of the analysis.

## Setup
### 1. Clone this repository 
### 2. Build conda environment
Change directories to the repository (`cd gfa_manuscript`). Use
```
conda env create -f gfa_manuscript_env.yml
```
to build the conda environment. This may take some time. I recommend doing this on a compute node with 10G of memory.

Type 
```
conda env list
```
to verify success. You should now see `gfa_manuscript` in your list of environments. 
Use 
```
conda activate gfa_manuscript
```
to activate the envrionment. 

## 3. Set up R package envrionments
Change directories to `code/applied_analysis`. Start R. 
In R use
```
renv::restore()
```
to install all necessary packages. This will take some time. 

Change directories to `code/simulations` and repeat this procedure. This
will go more quickly because you will not need to install as many new packages. 


## Reproduce the Applied analysis

### 1. Downloa data 

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
If you need to specify an account, you can do this on the --default-resources line, e.g. you 
may change this to

```
--default-resources mem_mb=5000 runtime=120 account=myaccount \
```

If your cluster does not use slurm, you will need to use a different executor. See 
the Snakemake documentation for more details. 

### 3. Execute the snakemake pipeline 

Use 

```
nohup ./run-snakemake-gfa_v9.sh &
```

Using nohup and the "&" will allow it to run uninterrupted as the process may take a few hours. 
Most  of the time is spent processing data and plotting. 

### 4. Look at reults. 
Check the end of nohup.out to make sure the pipeline ran without errors. If there was an error, 
check the listed log file to see why. You may need to change your specification of cluster resources. 
If the the snakemake job dies for some reason before completing, you can restart it from where it left off
by repeating the command in Step 3. At the
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
```
results/applied_analysis/gfa_results/metab_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS
``` 
and 
```
results/applied_analysis/gfa_results/bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS
``` 
All of the other results files were used for various supplementary comparisons. 


## Reproduce the simulations 

### 1. Change to the `code/simulations` directory. 
### 2. Modify run-snakemake-sims_v9.sh to work with your cluster. 
### 3. Execute the snakemake pipeline 

```
nohup ./run-snakemake-gfa_v9.sh &
```
