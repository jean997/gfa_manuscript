
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
    input = list(),
    output = list('../../results/simulation_results/simulate/simulate.B1.0.0.RDS', "out" = '../../results/simulation_results/simulate/simulate.B1.0.0.RDS'),
    params = list(),
    wildcards = list('B1', '0', '0', "scenario" = 'B1', "ndense" = '0', "rep" = '0'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list(),
    rule = 'simulate',
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
library(GWASBrewer)

source("R/rand_F_sim.R")


nblocks <- 10
block_size <- 5
ntrait <- 50
nfactor <- 12
nvar <- 5e5
N <- 50000
pi_L <- 1000/nvar
pi_theta <- 500/nvar
r2_thresh <- 0.01
sparse_mean <- 5

scenario <- snakemake@wildcards["scenario"]
ndense <- snakemake@wildcards["ndense"]
rep <- snakemake@wildcards["rep"]

set.seed(rep)

if(scenario == "A"){
    overlap_prop <- 0
    R <- diag(ntrait)
}else{
    overlap_prop = 1
    if(scenario == "B1"){
        R <- diag(ntrait)
    }else if(scenario == "B2"){
        Rblock <- matrix(0.3, block_size, block_size)
        diag(Rblock) <- 1
        R <- kronecker(diag(nblocks), Rblock);
    }else if(scenario == "B3"){
        Rblock <- matrix(0.7, block_size, block_size)
        diag(Rblock) <- 1
        R <- kronecker(diag(nblocks), Rblock);
    }
}

dat <- rand_F_sim(R, overlap_prop,
                  ntrait, nfactor,
                  ndense, sparse_mean,
                  N, nvar, pi_L, pi_theta);
#dat$ld_rand <- readRDS("ld_rand_5e5.RDS");
pmin <- apply(dat$pval, 1, min);
rmin <- apply(dat$rank, 1, min);
dat$ld_pval <- sim_ld_prune(dat, pvalue = pmin,
                            R_LD = ld_mat_list,
                            r2_thresh = r2_thresh,
                            pval_thresh = 1);

#dat$ld_rank <- sim_ld_prune(dat, pvalue = rmin,
#                            R_LD = ld_mat_list,
#                            r2_thresh = r2_thresh,
#                            pval_thresh = Inf);

saveRDS(dat, file = snakemake@output["out"])
