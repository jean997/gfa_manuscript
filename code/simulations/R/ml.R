library(nFactors)

dat <- readRDS(snakemake@input[["inp"]])
R <- readRDS(snakemake@input[["R"]])
cond_num <- as.numeric(snakemake@wildcards[["cond_num"]])

#nfact <- as.numeric(snakemake@wildcards[["nfct"]])
Rg <- GFA::condition(R$R_ldsc$Rg, cond_num = cond_num, corr = TRUE)

eR <- eigen(Rg, only.values=T)$values
nF <- nFactors::nScree(x = eR)$Components |> as.list()
i <- which(names(nF) == "nparallel")
nF <- nF[-i]
nF$six <- 6
nF$twelve <- 12

fits <- lapply(nF, function(nfact){
    if(nfact == 0){
        fit <- list("F_hat" = matrix(NA, nrow = ncol(Z_hat), ncol = 0), 
                    "L_hat" = NULL,
                     "ix" = ix) 
        return(fit)
    }

    f <- try(stats::factanal(covmat = Rg, factors = nfact, rotation="promax"), silent = TRUE);
    if(inherits(f, "try-error")){
        fit <- list("F_hat" = NULL,
                  "L_hat" = NULL, 
                  "B_hat" = NULL, 
                  "fit" = f);
    }else{
        fit <- list("F_hat" = matrix(f$loadings, nrow = nrow(Rg), ncol = nfact), 
                  "L_hat" = NULL, 
                  "fit" = f);
    }
    fit
    })

saveRDS(fits, file = snakemake@output[["out"]])

