library(GWASBrewer)

source("R/rand_F_sim.R")

scenario <- snakemake@wildcards[["scenario"]]
rep <- as.numeric(snakemake@wildcards[["rep"]])

set.seed(rep)

if(! scenario %in% c("bc", "bcshuffle")){
    ndense <- as.numeric(snakemake@wildcards[["ndense"]])
    nblocks <- 10
    block_size <- 5
    ntrait <- 50
    nfactor <- 12
    nvar <- 1e6
    N <- 50000
    pi_L <- 1000/nvar
    pi_theta <- 500/nvar
    r2_thresh <- 0.01
    sparse_mean <- 5

    if(scenario %in% paste0(c("A", "B1", "B2", "B3"), "null")){
        nfactor <- 0
        scenario <- stringr::str_replace(scenario, "null", "")
    }
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
    if(nfactor > 0){
        dat <- rand_F_sim(R, overlap_prop,
                  ntrait, nfactor,
                  ndense, sparse_mean,
                  N, nvar, pi_L, pi_theta);
    }else{
        h2_trait <- runif(n = ntrait, 0.05, 0.2)
        dat <- sim_mv(G = ntrait, N = N*matrix(1, nrow = ntrait, ncol = ntrait), 
                      J = nvar, 
                  h2 = h2_trait,
                  pi = pi_L + pi_theta, 
                  est_s = TRUE, 
                  R_E = R,
                  af = AF, 
                  R_LD = ld_mat_list)
        dat$pval <- with(dat, 2*pnorm(-abs(beta_hat/s_estimate)))
        dat$rank <- apply(dat$pval, 2, rank)
        dat$N <- rep(N, ntrait)
        dat <- GWASBrewer:::calc_ld_scores(dat, R_LD = ld_mat_list)
    }
}else{
    R <- readRDS("bc_R_estimate.R_ldsc.RDS")
    RE <- readRDS("bc_RE.RDS")
    ntrait <- ncol(RE)
    h2_trait <- diag(R$Rg)

    N <- 166000
    overlap_prop <- 1
    myF <- readRDS("bc_Fhat.RDS")


    nfactor <- ncol(myF)
    if(scenario == "bcshuffle"){
        myF <- myF[sample(seq(ntrait), size = ntrait, replace = FALSE),] 
    }
    ## in estimate from real bloodcell data omega = rowSums(pve) is 1 or very close for all traits
    omega <- rep(1, ntrait)
    ## we will set h2_factor to 1 since we are inputting the true R_obs
    h2_factor <- rep(1, nfactor)
    myN <- matrix(N, nrow = ntrait, ncol = ntrait)

    nvar <- 1e6
    #pi_L <- rep(1000/nvar, nfactor)
    #pi_theta <- 500/nvar # not used because omega is 1
    r2_thresh <- 0.01

    snp_eff_fcts <- readRDS("bc_factor_dists.RDS")

    dat <- sim_lf(F_mat = myF,
                  N=myN,
                  J=nvar,
                  h2_trait = h2_trait,
                  omega = omega,
                  h2_factor = h2_factor,
                  pi_L = 1,
                  pi_theta = 1,
                  R_E = RE,
                  R_LD = ld_mat_list,
                  af = AF,
                  snp_effect_function_L = snp_eff_fcts,
                  est_s = TRUE)
     dat$pval <- with(dat, 2*pnorm(-abs(beta_hat/s_estimate)))
     dat$rank <- apply(dat$pval, 2, rank)
     dat$N <- rep(N, ntrait)
     dat <- GWASBrewer:::calc_ld_scores(dat, R_LD = ld_mat_list)
}

#dat$nsig_per_trait <- sapply(1:ntrait, function(i){
#                                  x <- sim_ld_prune(dat, pvalue = i, 
#                                                    R_LD = ld_mat_list, 
#                                                    r2_thresh = 0.01, 
#                                                    pval_thresh = 5e-8)
#                                  cat(length(x), "\n")
#                                  return(length(x))})
#




dat$ld_rand <- readRDS("ld_rand_1e6.RDS");
pmin <- apply(dat$pval, 1, min);
#rmin <- apply(dat$rank, 1, min);
dat$ld_pval <- sim_ld_prune(dat, pvalue = pmin,
                            R_LD = ld_mat_list,
                            r2_thresh = r2_thresh,
                            pval_thresh = 1);

#dat$ld_rank <- sim_ld_prune(dat, pvalue = rmin,
#                            R_LD = ld_mat_list,
#                            r2_thresh = r2_thresh,
#                            pval_thresh = Inf);
dat$nsig_per_trait <- colSums(dat$pval[dat$ld_pval,] < 5e-8)

saveRDS(dat, file = snakemake@output[["out"]])
saveRDS(dat$F_mat, file = snakemake@output[["F_out"]])

