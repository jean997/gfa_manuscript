
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
    input = list('../../data/gfa_intermediate_data/bc_zmat.10.RDS', '../../results/applied_analysis/gfa_results/bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS', '../../data/gfa_intermediate_data/bc_R_estimate.R_ldsc.RDS', "z_file" = '../../data/gfa_intermediate_data/bc_zmat.10.RDS', "gfa_file" = '../../results/applied_analysis/gfa_results/bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS', "R" = '../../data/gfa_intermediate_data/bc_R_estimate.R_ldsc.RDS'),
    output = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS', "out" = '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS'),
    params = list(),
    wildcards = list('bc', '1_cn100', 'r20.01_kb1000_pvalue_jitter0_0', 'ldsc', '10', "prefix" = 'bc', "fs" = '1_cn100', "ldstring" = 'r20.01_kb1000_pvalue_jitter0_0', "Rstring" = 'ldsc', "chrom" = '10'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "jitter_sd" = c(0, 0.5), "jitter_seed" = c(1, 2, 3)), "R" = list("type" = c('ldsc', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = c(100, 1000)), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1), "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
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

gls_sol <- gfa_loadings_gls(Z_hat, S, gfafit)

gls_zscores <- gls_sol$L/gls_sol$S
gls_pvals <- gls_sol$P

nf <- ncol(gfafit$F_hat)


first_pos_name <- str_subset(names(data), "\\.pos$")[1]

res <- data.frame(cbind(gls_zscores, gls_pvals))
names(res) <- c(paste0("factor", 1:nf, ".z"), paste0("factor", 1:nf, ".p"))
res <- bind_cols(data[,c("chrom", "snp", first_pos_name,  "REF", "ALT")], res)
names(res)[3] <- "pos"
saveRDS(res, file = out)
