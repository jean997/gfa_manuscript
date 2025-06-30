library(GFA)
library(GWASBrewer)

fo <- readRDS("../../results/applied_analysis/gfa_results/bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS")

N <- as.numeric(fo$scale^2)
p <- length(N)
Nmat <- diag(N)
for(i in 1:p){
    for(j in 1:p){
        Nmat[i,j] <- min(N[i], N[j])
    }
}

R <- readRDS("../../data/gfa_intermediate_data/bc_R_estimate.R_ldsc.RDS")
h2_trait <- diag(R$Rg)

omega <- rowSums(fo$gfa_pve$pve[, -c(13:14)])
myF <- fo$F_hat_multi
nfactor <- ncol(myF)

Robs <- cov2cor(R$R)


nvar <- 5e5
pi_L <- 1000/nvar
pi_theta <- 500/nvar
r2_thresh <- 0.01
set.seed(1)
#h2_factor <- runif(n = ncol(myF), min = 0.1, max = 0.9)
h2_factor <- rep(1, nfactor)

dat <- sim_lf(F_mat = myF,
              N=Nmat,
              J=nvar,
              h2_trait = h2_trait,
              omega = omega,
              h2_factor = h2_factor,
              pi_L = rep(pi_L, nfactor),
              pi_theta = pi_theta,
              R_obs = Robs,
              R_LD = ld_mat_list,
              af = AF,
              est_s = TRUE)


s <- sample(1:p, size = p, replace = F)
myF <- fo$F_hat_multi[s,]
dat <- sim_lf(F_mat = myF,
              N=Nmat,
              J=nvar,
              h2_trait = h2_trait[s],
              omega = omega[s],
              h2_factor = h2_factor,
              pi_L = rep(pi_L, nfactor),
              pi_theta = pi_theta,
              R_E = Robs,
              R_LD = ld_mat_list,
              af = AF,
              est_s = TRUE)

dat$pval <- with(dat, 2*pnorm(-abs(beta_hat/s_estimate)))
pmin <- apply(dat$pval, 1, min);
dat$N <- N
dat <- GWASBrewer:::calc_ld_scores(dat, R_LD = ld_mat_list)

dat$ld_pval <- sim_ld_prune(dat, pvalue = pmin,
                            R_LD = ld_mat_list,
                            r2_thresh = r2_thresh,
                            pval_thresh = 1);

