source("renv/activate.R")
library(nFactors)

dat <- readRDS(snakemake@input[["inp"]])
pval_thresh <- as.numeric(snakemake@wildcards[["pthresh"]])
zero_thresh <- as.numeric(snakemake@wildcards[["zero_thresh"]])

pmin <- apply(dat$pval[dat$ld_pval,], 1, min);
ix <- dat$ld_pval[pmin < pval_thresh];

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,]);
Z_hat[abs(Z_hat) < zero_thresh] <- 0;

s <- svd(Z_hat);
pve <- with(s, d^2/sum(d^2));

nF <- nScree(x = s$d)
nF <- as.list(nF$Components)
i <- which(names(nF) == "nparallel")
nF <- nF[-i]
nF$six <- 6
nF$twelve <- 12

fit <- list("F_hat" = s$v, "L_hat" = with(s, u %*% diag(d)),
             "ix" = ix, 
             "n_factors" = nF)

########### Select factors
#ntrait <- ncol(dat$beta_hat)
#fct_ix <- which(fito$pve > 1/ntrait);
#u <- s$u[,fct_ix, drop = FALSE];
#v <- s$v[,fct_ix, drop = FALSE];
#d <- s$d[fct_ix];
#pve <- fito$pve[fct_ix];
#fit <- list("F_hat" = v, "L_hat" = u %*% diag(d), "B_hat" = u %*% diag(d) %*% t(v), "fit" = s, "ix" = fito$ix);


saveRDS(fit, file = snakemake@output[["out"]])
