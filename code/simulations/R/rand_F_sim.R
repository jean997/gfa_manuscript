library(GWASBrewer);

rand_F_sim <- function(R, overlap_prop,
                       ntrait, nfactor,
                       ndense, sparse_mean,
                       N, nvar, pi_L, pi_theta){

  g_F <- function(n){runif(n, -1, 1)}

  nz_factor <- c(pmin(rpois(ndense, (ntrait/2)-1)+1, ntrait), pmin( rpois(nfactor-ndense, sparse_mean - 2)+2, ntrait))

  h2_trait <- runif(n=ntrait, 0.05, 0.20)
  omega <- runif(n=ntrait, 0.2, 0.7)
  h2_factor <- runif(n=nfactor, 0.3, 0.8)


  data("ld_mat_list")
  data("AF")
  done <- FALSE
  while(!done){
    myF <- generate_random_F(K = nfactor, M = ntrait,
                             g_F = g_F, nz_factor = nz_factor,
                             omega = omega, h2_trait = h2_trait,
                             pad = FALSE)
    if(overlap_prop == 0){
      myN = N
    }else if(overlap_prop == 1){
      myN = matrix(N, nrow = ntrait, ncol = ntrait)
    }

    missing_traits <- which(rowSums(myF^2) == 0)
    if(length(missing_traits) > 0){
      omega[missing_traits] <- 0
    }
    dat <- try(sim_lf(F_mat = myF,
                      N=myN,
                      J=nvar,
                      h2_trait = h2_trait,
                      omega = omega,
                      h2_factor = h2_factor,
                      pi_L = rep(pi_L, nfactor),
                      pi_theta = pi_theta,
                      R_E = R,
                      R_LD = ld_mat_list,
                      af = AF,
                      est_s = TRUE), silent = TRUE)
    if(!inherits(dat, "try-error")){ done <- TRUE}
  }
  dat$pval <- with(dat, 2*pnorm(-abs(beta_hat/s_estimate)))
  dat$rank <- apply(dat$pval, 2, rank)
  dat$N <- rep(N, ntrait)
  dat <- GWASBrewer:::calc_ld_scores(dat, R_LD = ld_mat_list)
  return(dat)
}
