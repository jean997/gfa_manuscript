
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
    input = list('../../data/gfa_intermediate_data/bc90_zmat.1.RDS', '../../data/gfa_intermediate_data/bc90_zmat.2.RDS', '../../data/gfa_intermediate_data/bc90_zmat.3.RDS', '../../data/gfa_intermediate_data/bc90_zmat.4.RDS', '../../data/gfa_intermediate_data/bc90_zmat.5.RDS', '../../data/gfa_intermediate_data/bc90_zmat.6.RDS', '../../data/gfa_intermediate_data/bc90_zmat.7.RDS', '../../data/gfa_intermediate_data/bc90_zmat.8.RDS', '../../data/gfa_intermediate_data/bc90_zmat.9.RDS', '../../data/gfa_intermediate_data/bc90_zmat.10.RDS', '../../data/gfa_intermediate_data/bc90_zmat.11.RDS', '../../data/gfa_intermediate_data/bc90_zmat.12.RDS', '../../data/gfa_intermediate_data/bc90_zmat.13.RDS', '../../data/gfa_intermediate_data/bc90_zmat.14.RDS', '../../data/gfa_intermediate_data/bc90_zmat.15.RDS', '../../data/gfa_intermediate_data/bc90_zmat.16.RDS', '../../data/gfa_intermediate_data/bc90_zmat.17.RDS', '../../data/gfa_intermediate_data/bc90_zmat.18.RDS', '../../data/gfa_intermediate_data/bc90_zmat.19.RDS', '../../data/gfa_intermediate_data/bc90_zmat.20.RDS', '../../data/gfa_intermediate_data/bc90_zmat.21.RDS', '../../data/gfa_intermediate_data/bc90_zmat.22.RDS', 'orig_bloodcell_ldsc90.csv', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.ldscore.gz', "Z" = c('../../data/gfa_intermediate_data/bc90_zmat.1.RDS', '../../data/gfa_intermediate_data/bc90_zmat.2.RDS', '../../data/gfa_intermediate_data/bc90_zmat.3.RDS', '../../data/gfa_intermediate_data/bc90_zmat.4.RDS', '../../data/gfa_intermediate_data/bc90_zmat.5.RDS', '../../data/gfa_intermediate_data/bc90_zmat.6.RDS', '../../data/gfa_intermediate_data/bc90_zmat.7.RDS', '../../data/gfa_intermediate_data/bc90_zmat.8.RDS', '../../data/gfa_intermediate_data/bc90_zmat.9.RDS', '../../data/gfa_intermediate_data/bc90_zmat.10.RDS', '../../data/gfa_intermediate_data/bc90_zmat.11.RDS', '../../data/gfa_intermediate_data/bc90_zmat.12.RDS', '../../data/gfa_intermediate_data/bc90_zmat.13.RDS', '../../data/gfa_intermediate_data/bc90_zmat.14.RDS', '../../data/gfa_intermediate_data/bc90_zmat.15.RDS', '../../data/gfa_intermediate_data/bc90_zmat.16.RDS', '../../data/gfa_intermediate_data/bc90_zmat.17.RDS', '../../data/gfa_intermediate_data/bc90_zmat.18.RDS', '../../data/gfa_intermediate_data/bc90_zmat.19.RDS', '../../data/gfa_intermediate_data/bc90_zmat.20.RDS', '../../data/gfa_intermediate_data/bc90_zmat.21.RDS', '../../data/gfa_intermediate_data/bc90_zmat.22.RDS'), "gwas_info" = 'orig_bloodcell_ldsc90.csv', "m" = c('/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.M_5_50'), "l2" = c('/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.ldscore.gz')),
    output = list('../../data/gfa_intermediate_data/bc90_R_estimate.R_ldsc.RDS', "out" = '../../data/gfa_intermediate_data/bc90_R_estimate.R_ldsc.RDS'),
    params = list('1e3', "cond_num" = '1e3'),
    wildcards = list('bc90', "prefix" = 'bc90'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 4758, "mem_mib" = 4538, "disk_mb" = 4758, "disk_mib" = 4538, "tmpdir" = '/tmp'),
    config = list("input" = list("bc95" = 'orig_bloodcell_ldsc95.csv', "bc97" = 'orig_bloodcell_ldsc97.csv', "bc90" = 'orig_bloodcell_ldsc90.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = '1e3'), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1, 2, 3), "af_thresh" = 0.01, "sample_size_tol" = 0.1, "cor_clust" = 1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'R_ldsc_full',
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
library(readr)
library(GFA)
library(stringr)


# sample size affects genetic covariance and h2 but not intercept or genetic correlation
gwas_info <- read_csv(snakemake@input[["gwas_info"]])
z_files = unlist(snakemake@input[["Z"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]
cond_num = as.numeric(snakemake@params[["cond_num"]])

if(!is.finite(cond_num)){
  mwc <- FALSE
}else{
  mwc <- TRUE
}

ld <- purrr::map_dfr(1:22, function(c){
  read_table(ld_files[c])
})

M <- purrr:::map(1:22, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()

X <- map_dfr(z_files, function(f){
  readRDS(f) %>%
    rename(SNP = snp) %>%
    inner_join(., ld)})

names <- gwas_info$name

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()
nmsz <- str_replace(colnames(Z_hat), ".z$", "")
SS <- X %>%
  select(ends_with(".ss")) %>%
  as.matrix()
nmss <- str_replace(colnames(SS), ".ss$", "")
o <- match(nmsz, nmss)
SS <- SS[, o]
N <- apply(SS, 2, median)
if(any(is.na(N))){
  N[is.na(N)] <- gwas_info$pub_sample_size[is.na(N)]
}

R <- R_ldsc(Z_hat = Z_hat,
              ldscores = X$L2,
              ld_size = M,
              N = N,
              return_gencov = TRUE,
              make_well_conditioned = mwc,
              cond_num = cond_num
              )

ret <- list(R = R$Se, Rg = R$Sg, names = nmsz)
saveRDS(ret, file=out)

