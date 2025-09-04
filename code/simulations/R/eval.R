source("renv/activate.R")
library(dplyr)
library(GFA)


eval <- function(fit, f_true, n_subset = NULL){
    disc_thresh = c(0.9, 0.95, 0.98);
    extra_thresh = c(0.8, 0.7);


    if(is.null(fit$F_hat)){
        disc <- NULL;
        n_disc <- rep(NA, length(disc_thresh));
        n_extra <- rep(NA, length(extra_thresh));
        frob_n <- NA 
        n_total <- NA 
    }else{
        if(!inherits(fit$F_hat, "matrix")){
            fit$F_hat <- matrix(fit$F_hat, ncol=1)
        }
        nfactors <- ncol(fit$F_hat)
        if(is.null(n_subset)){
            n_subset <- nfactors
        }
        if(n_subset == 0 | nfactors == 0){
            n_disc <- rep(0, length(disc_thresh))
            n_extra <- rep(0, length(extra_thresh))
            n_total <- 0
            if(is.null(f_true)){
                   ntrue <- 0
            }else{
                ntrue <- ncol(f_true) - length(GFA:::find_single_trait(f_true))
            }
            frob_n <- sqrt(ntrue)
            res <- data.frame(frob_n = frob_n, 
                      n_extra_0.8 = n_extra[1],
                      n_extra_0.7 = n_extra[2],
                      n_disc_0.9  = n_disc[1],                 
                      n_disc_0.95  = n_disc[2],
                      n_disc_0.98   = n_disc[3], 
                      n_total = n_total)
            return(res)
         }
         f_hat <- fit$F_hat[, 1:n_subset, drop = FALSE]
         if(is.null(f_true)){
            f_true <- diag(nrow(fit$F_hat))
         }
         sol <- min_norm(f_true = f_true, f_hat = f_hat, single_trait_thresh = 0.98)
        
        n_disc <- sapply(disc_thresh, function(x){sum(sol$solution$match_score > x, na.rm=T)});
        n_extra <- sapply(extra_thresh, function(x){sum(sol$solution$match_score < x & !is.na(sol$solution$est_ix), na.rm=T)});
        n_total <- sum(!is.na(sol$solution$est_ix))
        frob_n <- sol$frob_n
    }
    res <- data.frame(frob_n = frob_n, 
                      n_extra_0.8 = n_extra[1],
                      n_extra_0.7 = n_extra[2],
                      n_disc_0.9  = n_disc[1],                 
                      n_disc_0.95  = n_disc[2],
                      n_disc_0.98   = n_disc[3], 
                      n_total = n_total)
    return(res)
}

f_true <- readRDS(snakemake@input[["f_true"]])
fit <- readRDS(snakemake@input[["fit"]])



if(snakemake@wildcards[["method"]]  %in%  c("ml", "guide", "factorgo")){
   df <- purrr::map_dfr(seq(length(fit)), function(i){
                dff <- eval(fit[[i]], f_true)
                dff$subset <- names(fit)[i]
                dff
                      })
}else if(is.null(fit$n_factors)){
    df <- eval(fit, f_true)
    df$subset <- "all"
}else{
    df <- purrr::map_dfr(1:length(fit$n_factors), function(i){
                             dff <- eval(fit, f_true, n_subset = fit$n_factors[[i]])
                             dff$subset <- names(fit$n_factors)[i]
                             dff
                      })
}

print(df)

df$method <- paste0(snakemake@wildcards[["method"]], "-", snakemake@wildcards[["mtype"]])
df$scenario <- snakemake@wildcards[["scenario"]]
df$ndense <- snakemake@wildcards[["ndense"]]
df$rep <- snakemake@wildcards[["rep"]]
df$ftrue_file <- snakemake@input[["f_true"]]
df$fit_file <- snakemake@input[["fit"]]

#obj <- list(df = df, min_norm_func = GFA::min_norm, eval_func = eval, f_true = f_true, fit = fit)
saveRDS(df, file = snakemake@output[["out"]])




