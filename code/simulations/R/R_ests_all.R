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
