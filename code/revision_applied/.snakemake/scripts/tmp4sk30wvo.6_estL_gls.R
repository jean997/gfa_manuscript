
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
    input = list('../../data/gfa_intermediate_data/metab_zmat.3.RDS', '../../results/gfa_results/metab_gfa_gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_pt0.05.final.RDS', '../../data/gfa_intermediate_data/metab_R_estimate.ldpruned_r20.01_kb1000_pvalue.R_pt0.05.RDS', "z_file" = '../../data/gfa_intermediate_data/metab_zmat.3.RDS', "gfa_file" = '../../results/gfa_results/metab_gfa_gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_pt0.05.final.RDS', "R" = '../../data/gfa_intermediate_data/metab_R_estimate.ldpruned_r20.01_kb1000_pvalue.R_pt0.05.RDS'),
    output = list('../../results/gfa_results/metab_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_pt0.05.3.RDS', "out" = '../../results/gfa_results/metab_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_pt0.05.3.RDS'),
    params = list(),
    wildcards = list('metab', '1', 'r20.01_kb1000_pvalue', 'pt0.05', '3', "prefix" = 'metab', "fs" = '1', "ldstring" = 'r20.01_kb1000_pvalue', "Rstring" = 'pt0.05', "chrom" = '3'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = 100), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = 1, "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/gfa_results/')),
    rule = 'gls_loadings_chrom',
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
data <- readRDS(snakemake@input[["z_file"]])
gfafit <- readRDS(snakemake@input[["gfa_file"]])
R <- readRDS(snakemake@input[["R"]])
stopifnot(all(R$names == gfafit$names))
R <- R$R

Z_hat <- data %>%
  select(ends_with(".z")) %>%
  as.matrix()

dat_names <- str_replace(colnames(Z_hat), ".z$", "")
ix <- which(dat_names %in% gfafit$names)
Z_hat <- Z_hat[,ix]
dat_names <- dat_names[ix]

o <- match(dat_names, gfafit$names)
Z_hat <- Z_hat[,o]
S <- matrix(1, nrow = nrow(Z_hat), ncol = ncol(Z_hat))
gls_sol <- GFA:::loadings_gls(Z_hat, S, R, gfafit$F_hat)
gls_zscores <- gls_sol$L/gls_sol$S
gls_pvals <- 2*pnorm(-abs(gls_zscores))
nf <- ncol(gfafit$F_hat)
res <- data.frame(cbind(gls_zscores, gls_pvals))
names(res) <- c(paste0("factor", 1:nf, ".z"), paste0("factor", 1:nf, ".p"))
res <- bind_cols(data[,c("chrom", "snp", "REF", "ALT")], res)
saveRDS(res, file = out)
