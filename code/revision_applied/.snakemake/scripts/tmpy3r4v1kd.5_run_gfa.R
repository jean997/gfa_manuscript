
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
    input = list('../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.1.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.2.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.3.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.4.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.5.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.6.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.7.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.8.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.9.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.10.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.11.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.12.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.13.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.14.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.15.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.16.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.17.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.18.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.19.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.20.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.21.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.22.RDS', '../../data/gfa_intermediate_data/bc97v2_R_estimate.ldpruned_r20.01_kb1000_pvalue.R_pt0.05_cc1.RDS', "Z" = c('../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.1.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.2.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.3.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.4.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.5.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.6.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.7.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.8.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.9.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.10.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.11.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.12.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.13.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.14.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.15.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.16.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.17.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.18.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.19.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.20.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.21.RDS', '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.22.RDS'), "R" = '../../data/gfa_intermediate_data/bc97v2_R_estimate.ldpruned_r20.01_kb1000_pvalue.R_pt0.05_cc1.RDS'),
    output = list('../../results/applied_analysis/gfa_results/bc97v2_gfa_gfaseed3.ldpruned_r20.01_kb1000_pvalue.R_pt0.05_cc1.1.RDS', "out" = '../../results/applied_analysis/gfa_results/bc97v2_gfa_gfaseed3.ldpruned_r20.01_kb1000_pvalue.R_pt0.05_cc1.1.RDS'),
    params = list('default', 'Inf', "params_file" = 'default', "max_snps" = 'Inf'),
    wildcards = list('bc97v2', '3', 'r20.01_kb1000_pvalue', 'pt0.05_cc1', "prefix" = 'bc97v2', "fs" = '3', "ldstring" = 'r20.01_kb1000_pvalue', "Rstring" = 'pt0.05_cc1'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("input" = list("bc95v2" = 'orig_bloodcell_ldsc95.csv', "bc97v2" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = 100), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1, 2, 3), "af_thresh" = 0.01, "sample_size_tol" = 0.1, "cor_clust" = 1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
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

z_files = unlist(snakemake@input[["Z"]])



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

# if(str_ends(R_est_file, "none_R.txt")){
#   R <- list(names = nms, R = diag(length(nms), nrow = ntrait))
# }else{
R <- readRDS(R_est_file)
stopifnot(all(R$names %in% nms))
z_order <- match(R$names, nms)
SS <- SS[,z_order]
Z_hat <- Z_hat[,z_order]
#R$R <- cov2cor(R$R)
#}

rownames(R$R) <- colnames(R$R) <- NULL


N <- apply(SS, 2, median)
t <- system.time(f <- gfa_fit(Z_hat = Z_hat,
                                N = N,
                                R = R$R,
                                params = params,
                                mode = "z-score",
                                method = "fixed_factors"))


f$snps <- snps
f$names <- R$names
f$time <- t
saveRDS(f, file=out)

