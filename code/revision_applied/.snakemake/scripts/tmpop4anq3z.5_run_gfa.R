
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
    input = list('../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.1.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.2.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.3.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.4.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.5.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.6.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.7.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.8.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.9.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.10.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.11.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.12.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.13.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.14.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.15.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.16.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.17.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.18.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.19.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.20.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.21.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.22.RDS', '../../data/gfa_intermediate_data/bc_R_estimate.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_pt0.05.RDS', 'orig_bloodcell_ldsc97.csv', "Z" = c('../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.1.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.2.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.3.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.4.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.5.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.6.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.7.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.8.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.9.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.10.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.11.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.12.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.13.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.14.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.15.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.16.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.17.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.18.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.19.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.20.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.21.RDS', '../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.22.RDS'), "R" = '../../data/gfa_intermediate_data/bc_R_estimate.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_pt0.05.RDS', "gwas_info" = 'orig_bloodcell_ldsc97.csv'),
    output = list('../../results/applied_analysis/gfa_results/bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_pt0.05.1.RDS', "out" = '../../results/applied_analysis/gfa_results/bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_pt0.05.1.RDS'),
    params = list('default', 'Inf', "params_file" = 'default', "max_snps" = 'Inf'),
    wildcards = list('bc', '1', '100', 'r20.01_kb1000_pvalue_jitter0.5_1', 'pt0.05', "prefix" = 'bc', "fs" = '1', "cond_num" = '100', "ldstring" = 'r20.01_kb1000_pvalue_jitter0.5_1', "Rstring" = 'pt0.05'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "jitter_sd" = c(0, 0.5), "jitter_seed" = c(1, 2, 3)), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = c(100)), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1, 2, 3), "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'run_gfa',
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


out <- snakemake@output[["out"]]
R_est_file <- snakemake@input[["R"]]
params_file <- snakemake@params[["params_file"]]
max_snp <- as.numeric(snakemake@params[["max_snps"]])
seed <- snakemake@wildcards[["fs"]]
info <- readr::read_csv(snakemake@input[["gwas_info"]])

z_files = unlist(snakemake@input[["Z"]])

cond_num <- as.numeric(snakemake@wildcards[["cond_num"]])


set.seed(seed)



if(params_file == "default"){
  params <- gfa_default_parameters()
}else{
  params <- readRDS(params_file)
}


# Read in data
X <- map_dfr(z_files, readRDS)

ntrait <- X %>%
          select(ends_with(".z")) %>%
          ncol()


if(nrow(X) > max_snp){
    ix <- sample(seq(nrow(X)), size = max_snp, replace = FALSE)
    X <- X[ix,]
}

Z_hat <- X %>%
         select(ends_with(".z")) %>%
         as.matrix()

SS <- X %>%
      select(ends_with(".ss")) %>%
      as.matrix()

snps <- X$snp

nms <- names(X)[grep(".z$", names(X))] %>% str_replace(".z$", "")

R <- readRDS(R_est_file)


## condition
#d <- diag(R$R)
#corR <- cov2cor(R$R)
#newR <- GFA::condition(R = corR,  cond_num=cond_num, corr = TRUE)
#newR <- diag(sqrt(d)) %*% newR %*% diag(sqrt(d))
newR <- Matrix::nearPD(R$R, posd.tol = 1/cond_num, keepDiag = TRUE)
newR <- newR$mat


stopifnot(all(R$names %in% nms))
z_order <- match(R$names, nms)
SS <- SS[,z_order]
Z_hat <- Z_hat[,z_order]

rownames(newR) <- colnames(newR) <- NULL


N <- apply(SS, 2, median)

## get case count and population prevalence
if(!"case_count" %in% names(info)){
    N_case = NULL
}else if(all(is.na(info$case_count))){
    N_case = NULL
}else{
    N_case <- info[match(R$names, info$name),]$case_count
}

if(!"pop_prev" %in% names(info)){
    pop_prev = NULL
}else if(all(is.na(info$pop_prev))){
    pop_prev = NULL
}else{
    pop_prev <- info[match(R$names, info$name),]$pop_prev
}


t <- system.time(f <- gfa_fit(Z_hat = Z_hat,
                                N = N,
                                N_case = N_case,
                                pop_prev = pop_prev,
                                R = newR, 
                                params = params,
                                mode = "z-score",
                                method = "fixed_factors"))


f$snps <- snps
f$names <- R$names
f$time <- t
saveRDS(f, file=out)

