library(GFA);

dat <- readRDS(snakemake@input[["inp"]])
Rtype <- snakemake@wildcards[["Rstring"]]
R <- readRDS(snakemake@input[["R"]])
cond_num <- as.numeric(snakemake@wildcards[["cond_num"]])

if(Rtype == "ldsc"){
    R <- R$R_ldsc$Se
}else if(Rtype == "pt0.05"){
    R <- R$R_pt
}else if(Rtype == "none"){
    R <- NULL
}else if(Rtype == "oracle"){
    R <- dat$R
}

#d <- diag(R)
#newR <- cov2cor(R)
#newR <- condition(R, cond_num = cond_num, corr = T)
#newR <- diag(sqrt(d)) %*% newR %*% diag(sqrt(d))
if(is.null(R)){
    newR <- R
}else{
    newR <- Matrix::nearPD(R, posd.tol = 1/cond_num, keepDiag = TRUE)
    newR <- as.matrix(newR$mat)
}


ldtype <- snakemake@wildcards[["ldtype"]]
pthresh <- as.numeric(snakemake@wildcards[["pthresh"]])

if(ldtype == "random"){
         ix <- readRDS("ld_rand_1e6.RDS")
}else if(ldtype == "pval"){
         ix <- dat$ld_pval;
}

pmin <- apply(dat$pval[ix,], 1, min);
ix <- ix[pmin <= pthresh]

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,]);
t <- system.time(
        fit <- gfa_fit(Z_hat = Z_hat, N = dat$N, 
                       R = newR));
fit$ix <- ix
fit$time <- t

saveRDS(fit, file = snakemake@output[["out"]])

