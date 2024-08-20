
library(GFA);
Z_hat <- with(dat, beta_hat/se_beta_hat);
ldscores <- dat$snp_info$l2;
myR <- GFA::R_ldsc(Z_hat, ldscores = ldscores, ld_size = nrow(Z_hat), 
                   N = dat$N, return_gencov=TRUE, 
                   make_well_conditioned = TRUE, cond_num = 1e3);
Rg <- GFA::condition(myR$Rg, 1e3, corr=TRUE);
R <- myR$Se;
