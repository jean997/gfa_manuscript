library(GFA);

dat <- readRDS(snakemake@input[["inp"]])
Rtype <- snakemake@wildcards[["Rstring"]]
R <- readRDS(snakemake@input[["R"]])
if(Rtype == "ldsc"){
    R <- R$R_ldsc$Se
}else if(Rtype == "pt0.05"){
    R <- R$R_pt
}else if(Rtype == "none"){
    R <- NULL
}else if(Rtype == "oracle"){
    R <- dat$R
}
ldtype <- snakemake@wildcards[["ldtype"]]
pthresh <- as.numeric(snakemake@wildcards[["pthresh"]])

if(ldtype == "random"){
         ix <- readRDS("ld_rand_5e5.RDS")
}else if(ldtype == "pval"){
         ix <- dat$ld_pval;
}

pmin <- apply(dat$pval[ix,], 1, min);
ix <- ix[pmin < pthresh]

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,]);
t <- system.time(
        fit <- gfa_fit(Z_hat = Z_hat, N = dat$N, 
                       R = R));
fit$ix <- ix;


saveRDS(fit, file = snakemake@output[["out"]])

