
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
    input = list('../../results/simulation_results/simulate/simulate.B1.3.0.RDS', "inp" = '../../results/simulation_results/simulate/simulate.B1.3.0.RDS'),
    output = list('../../results/simulation_results/R_ests/R_ests.B1.3.0.RDS', "out" = '../../results/simulation_results/R_ests/R_ests.B1.3.0.RDS'),
    params = list(),
    wildcards = list('B1', '3', '0', "scenario" = 'B1', "ndense" = '3', "rep" = '0'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1895, "mem_mib" = 1808, "disk_mb" = 1895, "disk_mib" = 1808, "tmpdir" = '/tmp'),
    config = list(),
    rule = 'estR',
    bench_iteration = as.numeric(NA),
    scriptdir = '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/R',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
library(GFA);

dat <- readRDS(snakemake@input[["inp"]])
## R-LDSC
Z_hat <- with(dat, beta_hat/se_beta_hat);
ldscores <- dat$snp_info$l2;
R_ldsc <- GFA::R_ldsc(Z_hat, ldscores = ldscores, ld_size = nrow(Z_hat), 
                   N = dat$N, return_gencov=TRUE, 
                   make_well_conditioned = TRUE, cond_num = 1e3);

R_ldsc$Rg <- cov2cor(R_ldsc$Rg)

### P-threshold
ix <- readRDS("ld_rand_5e5.RDS")
p_thresh <- 0.05
R_pt <- with(dat, R_pt(beta_hat[ix,], se_beta_hat[ix,], p_val_thresh = p_thresh, 
                                 make_well_conditioned = TRUE, cond_num = 1e3));

res <- list(R_ldsc = R_ldsc, R_pt = R_pt)
saveRDS(res, file = snakemake@output[["out"]])
