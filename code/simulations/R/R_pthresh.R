 R(library(GFA);
ix <- dat$ld_rand;
R <- with(dat, R_pt(beta_hat[ix,], se_beta_hat[ix,], p_val_thresh = p_thresh, 
                                 make_well_conditioned = TRUE, cond_num = 1e3));

