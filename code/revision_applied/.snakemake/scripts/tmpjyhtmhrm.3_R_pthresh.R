
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
    input = list('../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.1.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.2.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.3.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.4.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.5.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.6.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.7.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.8.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.9.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.10.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.11.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.12.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.13.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.14.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.15.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.16.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.17.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.18.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.19.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.20.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.21.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.22.RDS', "Z" = c('../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.1.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.2.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.3.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.4.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.5.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.6.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.7.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.8.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.9.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.10.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.11.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.12.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.13.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.14.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.15.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.16.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.17.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.18.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.19.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.20.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.21.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.22.RDS')),
    output = list('../../data/gfa_intermediate_data/bc_R_estimate.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_pt0.05.RDS', "out" = '../../data/gfa_intermediate_data/bc_R_estimate.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_pt0.05.RDS'),
    params = list(),
    wildcards = list('bc', 'r20.01_kb1000_pvalue_jitter0.5_1', '0.05', "prefix" = 'bc', "ldstring" = 'r20.01_kb1000_pvalue_jitter0.5_1', "pt" = '0.05'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "jitter_sd" = c(0, 0.5), "jitter_seed" = c(1, 2, 3)), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = c(100)), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1, 2, 3), "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'pt_R',
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
library(GFA)



p_thresh <- as.numeric(snakemake@wildcards[["pt"]])
out <- snakemake@output[["out"]]
z_files = unlist(snakemake@input[["Z"]])
#cond_num = as.numeric(snakemake@params[["cond_num"]])

#if(!is.finite(cond_num)){
#  mwc <- FALSE
#}else{
#  mwc <- TRUE
#}

# Read in data
X <- map_dfr(z_files, readRDS)

ntrait <- X %>%
  select(ends_with(".z")) %>%
  ncol()

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()

nms <- colnames(Z_hat) %>% stringr::str_replace(".z$", "")
Rpt <- R_pt(B_hat = Z_hat,
            S_hat = matrix(1, nrow = nrow(Z_hat), ncol = ncol(Z_hat)),
            p_val_thresh = p_thresh,
            return_cor = TRUE,
            make_well_conditioned = FALSE 
            )

ret <- list(R = Rpt, names = nms)
saveRDS(ret, file=out)
