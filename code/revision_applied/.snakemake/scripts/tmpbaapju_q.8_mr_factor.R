
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.22.RDS', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Bentham_2015_26502338/ebi-a-GCST003156.vcf.gz', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.Rgcor.RDS', "Z" = c('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.22.RDS'), "outcome_file" = '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Bentham_2015_26502338/ebi-a-GCST003156.vcf.gz', "R" = '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.Rgcor.RDS'),
    output = list('../../results/applied_analysis/mr_results/bc_mr_fct.Bentham_2015_26502338__lupus.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.RDS', "out" = '../../results/applied_analysis/mr_results/bc_mr_fct.Bentham_2015_26502338__lupus.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.RDS'),
    params = list('/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', 5e-08, 0.01, 1000, 'orig_autoimmune.csv', "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "pthresh" = 5e-08, "r2_thresh" = 0.01, "clump_kb" = 1000, "outcome_info_file" = 'orig_autoimmune.csv'),
    wildcards = list('bc', 'Bentham_2015_26502338__lupus', 'gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc', "prefix" = 'bc', "outcome" = 'Bentham_2015_26502338__lupus', "analysis" = 'gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 3354, "mem_mib" = 3199, "disk_mb" = 3354, "disk_mib" = 3199, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "jitter_sd" = c(0, 0.5), "jitter_seed" = c(1, 2, 3)), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = c(100)), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1, 2, 3), "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'mr_factor',
    bench_iteration = as.numeric(NA),
    scriptdir = '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_analysis/gfa_manuscript/code/revision_applied/R',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
library(ieugwasr)
library(MVMR)
library(dplyr)
library(ieugwasr)
library(purrr)
library(stringr)
library(gwasvcf)

source("R/format_ieu_chrom.R")

