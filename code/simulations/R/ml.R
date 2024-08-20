
dat <- readRDS(snakemake@input[["inp"]])
R <- readRDS(snakemake@input[["R"]])
nfact <- as.numeric(snakemake@wildcards[["nfct"]])
Rg <- R$R_ldsc$Rg

f <- try(stats::factanal(covmat = Rg, factors = nfact, rotation="promax"), silent = TRUE);
if(inherits(f, "try-error")){
        fit <- list("F_hat" = NULL,
                  "L_hat" = NULL, 
                  "B_hat" = NULL, 
                  "fit" = f);
}else{
        fit <- list("F_hat" = matrix(f$loadings, nrow = nrow(Rg), ncol = nfact), 
                  "L_hat" = NULL, 
                  "B_hat" = NULL, 
                  "fit" = f);
}

saveRDS(fit, file = snakemake@output[["out"]])

