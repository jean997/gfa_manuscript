
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
    input = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.22.RDS', "Z" = c('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.22.RDS')),
    output = list('../../results/applied_analysis/plots/bc_mhtplot.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.png', "mht_out" = '../../results/applied_analysis/plots/bc_mhtplot.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.png'),
    params = list(),
    wildcards = list('bc', 'gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc', "prefix" = 'bc', "analysis" = 'gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 3053, "mem_mib" = 2912, "disk_mb" = 3053, "disk_mib" = 2912, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "jitter_sd" = c(0, 0.5), "jitter_seed" = c(1, 2, 3)), "R" = list("type" = c('ldsc', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = c(100, 1000)), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1), "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'mht',
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
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)


# sample size affects genetic covariance and h2 but not intercept or genetic correlation
z_files = unlist(snakemake@input[["Z"]])
mht_out_prefix <- snakemake@output[["mht_out"]] |> str_replace(".1.png$", "")

X <- map_dfr(z_files, function(f){
  readRDS(f) 
})

X$chrom <- as.numeric(X$chrom)

fct_names <- str_subset(names(X), "\\.p$")
fct_names <- str_replace(fct_names, "\\.p$", "")

# order as shown in publication figures
if(str_detect(mht_out_prefix, "metab_")){
    fct_order <- c(1, 15, 5, 2, 10, 7,  4, 3, 13, 11,  6, 9,  8, 12, 14)
    disp_name <- paste0("Factor ", 1:14)
}else if(str_detect(mht_out_prefix, "bc_")){
    fct_order <- c(7, 2, 1,3,6, 14, 10, 9, 5,11, 4, 8, 12, 13)
    disp_name <- paste0("Factor ", 1:14)
}




nplot <- ceiling(length(fct_names)/2)


i <- 1
for(j in 1:nplot){
    png(paste0(mht_out_prefix, ".", j, ".png"), 
        height = 8, width = 15, units = "in", res = 300)
    par(mfrow = c(2, 1))
    nmax <- min(2, length(fct_names)-i -1)
    X$pval <- pmax(X[[paset0(fct_names[fct_order[i]], ".p")]], 1e-100)
    for(k in 1:nmax){
        qqman::manhattan(x = X, 
                         chr = "chrom", 
                         bp = "pos", 
                         p = "pval",
                         snp = "snp",
                         main = disp_name[i])
        i <- i + 1
    }
    dev.off()
}





