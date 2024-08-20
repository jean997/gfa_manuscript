
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
    input = list('../../results/simulation_results/simulate/simulate.A.3.0.RDS', '../../results/simulation_results/R_ests/R_ests.A.3.0.RDS', "inp" = '../../results/simulation_results/simulate/simulate.A.3.0.RDS', "R" = '../../results/simulation_results/R_ests/R_ests.A.3.0.RDS'),
    output = list('../../results/simulation_results/gfa/gfa_none_random_1.A.3.0.RDS', "out" = '../../results/simulation_results/gfa/gfa_none_random_1.A.3.0.RDS'),
    params = list(),
    wildcards = list('none', 'random', '1', 'A', '3', '0', "Rstring" = 'none', "ldtype" = 'random', "pthresh" = '1', "scenario" = 'A', "ndense" = '3', "rep" = '0'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1895, "mem_mib" = 1808, "disk_mb" = 1895, "disk_mib" = 1808, "tmpdir" = '/tmp'),
    config = list(),
    rule = 'gfa_R',
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
Rtype <- snakemake@wildcards[["Rstring"]]
R <- readRDS(snakemake@input[["R"]])
if(Rtype == "ldsc"){
    R <- R$R_ldsc$Se
}else if(Rtype == "pt0.05"){
    R <- R$R_pt
}else if(Rtype === "none"){
    R <- NULL
}else if(Rtype == "oracle"){
    R <- dat$R
}
ldtype <- snakemake@wildcards[["ldtype"]]
pthresh <- as.numeric(snakemake@wildcards[["pthresh"]])

if(ldtype == "random"){
         ix <- readRDS("ld_rand_5e5.RDS")
}else if(ldtype == "pval"){
         ix <- dat$ld_pval;
}

pmin <- apply(dat$pval[ix,], 1, min);
ix <- ix[pmin < pthresh]

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,]);
t <- system.time(
        fit <- gfa_fit(Z_hat = Z_hat, N = dat$N, 
                       R = R));
fit$ix <- ix;

df <- eval(fit)
df$method <- paste0("GFA-R", Rtype, "-", ldtype, "-",  pthresh)

res <- list(fit = fit, df)

saveRDS(res, file = snakemake@output[["out"]])

