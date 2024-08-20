
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
    input = list('../../results/simulation_results/simulate/simulate.B2.3.0.RDS', "inp" = '../../results/simulation_results/simulate/simulate.B2.3.0.RDS'),
    output = list('../../results/simulation_results/svd/svd_0_5e-8.B2.3.0.RDS', "out" = '../../results/simulation_results/svd/svd_0_5e-8.B2.3.0.RDS'),
    params = list(),
    wildcards = list('0', '5e-8', 'B2', '3', '0', "zero_thresh" = '0', "pthresh" = '5e-8', "scenario" = 'B2', "ndense" = '3', "rep" = '0'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1895, "mem_mib" = 1808, "disk_mb" = 1895, "disk_mib" = 1808, "tmpdir" = '/tmp'),
    config = list(),
    rule = 'svd',
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
source("R/eval.R")



dat <- readRDS(snakemake@input[["inp"]])
pval_thresh <- as.numeric(snakemake@wildcards[["pthresh"]])
zero_thresh <- as.numeric(snakemake@wildcards[["zero_thresh"]])

pmin <- apply(dat$pval[dat$ld_pval,], 1, min);
ix <- dat$ld_pval[pmin < pval_thresh];

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,]);
Z_hat[abs(Z_hat) < zero_thresh] <- 0;

s <- svd(Z_hat);
pve <- with(s, d^2/sum(d^2));
fito <- list("F_hat" = s$v, "L_hat" = with(s, u %*% diag(d)),
             "B_hat" = with(s, u %*% diag(d) %*% t(v)), "pve" = pve,
             "ix" = ix);

########### Select factors
ntrait <- ncol(dat$beta_hat)
fct_ix <- which(fito$pve > 1/ntrait);
u <- s$u[,fct_ix, drop = FALSE];
v <- s$v[,fct_ix, drop = FALSE];
d <- s$d[fct_ix];
pve <- fito$pve[fct_ix];
fit <- list("F_hat" = v, "L_hat" = u %*% diag(d), "B_hat" = u %*% diag(d) %*% t(v), "fit" = s, "ix" = fito$ix);


df <- eval(fit)
df$svd_zt <- zero_thresh
df$svd_pt <- pval_thresh
if(zero_thresh > 0){
    df$method <- "SVD-HT"
}else{
    df$method <- "SVD"
}
res <- list(fit = fit, eval = df)

saveRDS(res, file = snakemake@output[["out"]])
