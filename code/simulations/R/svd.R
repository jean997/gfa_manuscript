
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


saveRDS(fit, file = snakemake@output[["out"]])
