
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
    input = list('../../results/simulation_results/simulate/simulate.B1.12.0.RDS', "inp" = '../../results/simulation_results/simulate/simulate.B1.12.0.RDS'),
    output = list('../../results/simulation_results/spc/spc_5e-8.B1.12.0.RDS', "out" = '../../results/simulation_results/spc/spc_5e-8.B1.12.0.RDS'),
    params = list(),
    wildcards = list('5e-8', 'B1', '12', '0', "pthresh" = '5e-8', "scenario" = 'B1', "ndense" = '12', "rep" = '0'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1928, "mem_mib" = 1839, "disk_mb" = 1928, "disk_mib" = 1839, "tmpdir" = '/tmp'),
    config = list(),
    rule = 'spc',
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
library(PMA)
source("R/eval.R")



dat <- readRDS(snakemake@input[["inp"]])
pval_thresh <- as.numeric(snakemake@wildcards[["pthresh"]])

pmin <- apply(dat$pval[dat$ld_pval,], 1, min);
ix <- dat$ld_pval[pmin < pval_thresh];

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,]);
cv <- SPC.cv(x = Z_hat);
k <- ncol(Z_hat);
s <- try(SPC(x = Z_hat, sumabsv = cv$bestsumabsv, K = k), silent = TRUE);
while(class(s) == "try-error"){
    k  <- k -1;
    s <- try(SPC(x = Z_hat, sumabsv = cv$bestsumabsv, K = k), silent = TRUE);
};
pve <- with(s, c(prop.var.explained[1], diff(prop.var.explained)));
fito <- list("F_hat" = s$v, "L_hat" = with(s, u %*% diag(d)),
                   "B_hat" = with(s, u %*% diag(d) %*% t(v)), "fit"  = s, "pve" = pve,
                   "ix" = ix);

########### Select factors
fct_ix <- which(fito$pve > 1/ntrait);
u <- s$u[,fct_ix, drop = FALSE];
v <- s$v[,fct_ix, drop = FALSE];
d <- s$d[fct_ix];
pve <- fito$pve[fct_ix];
fit <- list("F_hat" = v, "L_hat" = u %*% diag(d), "B_hat" = u %*% diag(d) %*% t(v), "fit" = s, "ix" = fito$ix);

df <- eval(fit)
df$spc_pt <- pval_thresh
df$method <- "SPC"

res <- list(fit = fit, eval = df)

saveRDS(res, file = snakemake@output[["out"]])
