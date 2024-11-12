
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
    input = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.22.RDS', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.ldscore.gz', "Z" = c('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.22.RDS'), "m" = c('/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.M_5_50', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.M_5_50'), "l2" = c('/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/1.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/2.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/3.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/4.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/5.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/6.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/7.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/8.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/9.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/10.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/11.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/12.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/13.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/14.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/15.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/16.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/17.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/18.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/19.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/20.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/21.l2.ldscore.gz', '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/22.l2.ldscore.gz')),
    output = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.Rgcor.RDS', "out" = '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.Rgcor.RDS'),
    params = list(),
    wildcards = list('bc', 'gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none', "prefix" = 'bc', "analysis" = 'gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 2653, "mem_mib" = 2531, "disk_mb" = 2653, "disk_mib" = 2531, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = 100), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = 1, "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'R_ldsc_gls',
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
z_files = unlist(snakemake@input[["Z"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]

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

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()
nms <- str_replace(colnames(Z_hat), ".z$", "")

R <- R_ldsc(Z_hat = Z_hat,
            ldscores = X$L2,
            ld_size = M,
            N = rep(1, ncol(Z_hat)),
            return_gencov = TRUE,
            make_well_conditioned = FALSE # not needed we only need Rg
)

saveRDS(R, file=out)

