library(PMA)
library(nFactors)

dat <- readRDS(snakemake@input[["inp"]])
pval_thresh <- as.numeric(snakemake@wildcards[["pthresh"]])

pmin <- apply(dat$pval[dat$ld_pval,], 1, min);
ix <- dat$ld_pval[pmin < pval_thresh];

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,]);
cv <- SPC.cv(x = Z_hat);
k <- ncol(Z_hat);
s <- try(SPC(x = Z_hat, sumabsv = cv$bestsumabsv, K = k), silent = TRUE);
while(class(s) == "try-error"){
    k  <- k -1;
    s <- try(SPC(x = Z_hat, sumabsv = cv$bestsumabsv, K = k), silent = TRUE);
};
pve <- with(s, c(prop.var.explained[1], diff(prop.var.explained)));

nF <- nFactors::nScree(x = s$d)
nF <- as.list(nF$Components)
i <- which(names(nF) == "nparallel")
nF <- nF[-i]
nF$six <- 6
nF$twelve <- 12

fit <- list("F_hat" = s$v, "L_hat" = with(s, u %*% diag(d)),
              "fit"  = s, "pve" = pve,
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
