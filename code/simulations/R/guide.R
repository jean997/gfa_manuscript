library(nFactors)
library(fastICA)

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

seeds <- round(runif(n = length(nF), min = 1, max = 1e7))
fits <- lapply(seq(length(nF)), function(i){
    set.seed(seeds[i])
    n <- nF[[i]]
    if(n == 0){
        fit <- list("F_hat" = matrix(NA, nrow = ncol(Z_hat), ncol = 0), 
                    "L_hat" = NULL,
                     "ix" = ix) 
        return(fit)
    }
    ica_res <- fastICA(X = Z_hat, n.comp = n)
    W_LT <- s$v[, 1:n, drop = F] %*% ica_res$W
    W_XL <- s$u[,1:n, drop = F] %*% ica_res$W
    Sstar <- t(ica_res$W) %*% diag(s$d[1:n], nrow = n) %*% ica_res$W


    fit <- list("F_hat" = W_LT, "L_hat" = W_XL %*% Sstar,
                "W_XL" = W_XL, "ica_seed" = seeds[i],
                "ix" = ix) 
    fit
})
names(fits) <- names(nF)
saveRDS(fits, file = snakemake@output[["out"]])
