
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
    input = list('orig_bloodcell_ldsc97.csv', "gwas_info" = 'orig_bloodcell_ldsc97.csv'),
    output = list('../../data/gfa_intermediate_data/bc_R_estimate.R_none.RDS', "out" = '../../data/gfa_intermediate_data/bc_R_estimate.R_none.RDS'),
    params = list(),
    wildcards = list('bc', "prefix" = 'bc'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = 100), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = 1, "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/gfa_results/')),
    rule = 'none_R',
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


gwas_info <- snakemake@input[["gwas_info"]]
cat(gwas_info, "\n")
names <- readr::read_csv(gwas_info)$name
n <- length(names)
ret <- list(R = diag(n, nrow = n), names = names)
saveRDS(ret, file = snakemake@output[["out"]])
