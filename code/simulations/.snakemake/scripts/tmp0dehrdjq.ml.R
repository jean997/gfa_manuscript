
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
    input = list('../../results/simulation_results/simulate/simulate.A.12.0.RDS', '../../results/simulation_results/R_ests/R_ests.A.12.0.RDS', "inp" = '../../results/simulation_results/simulate/simulate.A.12.0.RDS', "R" = '../../results/simulation_results/R_ests/R_ests.A.12.0.RDS'),
    output = list('../../results/simulation_results/ml/ml_12.A.12.0.RDS', "out" = '../../results/simulation_results/ml/ml_12.A.12.0.RDS'),
    params = list(),
    wildcards = list('12', 'A', '12', '0', "nfct" = '12', "scenario" = 'A', "ndense" = '12', "rep" = '0'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1928, "mem_mib" = 1839, "disk_mb" = 1928, "disk_mib" = 1839, "tmpdir" = '/tmp'),
    config = list(),
    rule = 'ml',
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

dat <- readRDS(snakemake@input[["inp"]])
R <- readRDS(snakemake@input[["R"]])
nfact <- as.numeric(snakemake@wildcards[["nfct"]])
Rg <- R$R_ldsc$Rg

f <- try(stats::factanal(covmat = Rg, factors = nfact, rotation="promax"), silent = TRUE);
if(inherits(f, "try-error")){
        fit <- list("F_hat" = NULL,
                  "L_hat" = NULL, 
                  "B_hat" = NULL, 
                  "fit" = f);
}else{
        fit <- list("F_hat" = matrix(f$loadings, nrow = nrow(Rg), ncol = nfact), 
                  "L_hat" = NULL, 
                  "B_hat" = NULL, 
                  "fit" = f);
}
df <- eval(fit)
df$method <- paste0("ML-", nfact)
res <- list(fit = fit, df)

saveRDS(res, file = snakemake@output[["out"]])

