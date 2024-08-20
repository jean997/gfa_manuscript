
library(GFA)

eval <- function(fit, dat){
    disc_thresh = c(0.9, 0.95, 0.98);
    extra_thresh = c(0.8, 0.7);
    if(is.null(fit$F_hat)){
        disc <- NULL;
        n_disc <- rep(NA, length(disc_thresh));
        n_extra <- rep(NA, length(extra_thresh));
        frob_n <- frob_n_l <- opt_frob_n <- NA;
        err <- NA;
    }else{
        if(!"matrix" %in% class(fit$F_hat)){
            fit$F_hat <- matrix(fit$F_hat, ncol=1)
            if(!is.null(fit$L_hat)) fit$L_hat <- matrix(fit$L_hat, ncol=1);
        }
        if(is.null(fit$L_hat)){
            err <- NA;
            sol <- min_norm(f_true = dat$F_mat, f_hat = fit$F_hat);
            frob_n_l <- NA;
        }else{
            ix <- fit$ix;
            B <- dat$beta_marg[ix,] - dat$theta_marg[ix,];
            if(is.null(fit$scale)){
                B_hat <- (fit$L_hat %*% t(fit$F_hat))*dat$s_estimate[ix,];
            }else{
                B_hat <- (fit$L_hat %*% t(fit$F_hat*fit$scale))*dat$s_estimate[ix,];
            }
            err <- sqrt(sum((B_hat-B)^2)/sum(B^2));
            sol <- min_norm(f_true = dat$F_mat, f_hat = fit$F_hat,
                                l_true = dat$L_mat_marg[ix,],
                                l_hat = fit$L_hat);
            frob_n_l <- sol$frob_n_l
        };
        n_disc <- sapply(disc_thresh, function(x){sum(sol$solution$val > x, na.rm=T)});
        n_extra <- sapply(extra_thresh, function(x){sum(sol$solution$val < x & !is.na(sol$solution$est_ix), na.rm=T)});
        frob_n <- sol$frob_n
    }
    res <- data.frame(frob_n = frob_n, 
                      n_extra_0.8 = n_extra[1],
                      n_extra_0.7 = n_extra[2],
                      n_disc_0.9  = n_disc[1],                 
                      n_disc_0.95  = n_disc[2],
                      n_disc_0.98   = n_disc[3])
    return(res)
}

dat <- readRDS(snakemake@input[["dat"]])
fit <- readRDS(snakemake@input[["fit"]])
df <- eval(fit, dat)
df$method <- paste0(snakemake@wildcards[["method"]], "-", snakemake@wildcards[["mtype"]])
df$scenario <- snakemake@wildcards[["scenario"]]
df$ndense <- snakemake@wildcards[["ndense"]]
df$rep <- snakemake@wildcards[["rep"]]
df$data_file <- snakemake@input[["dat"]]
df$fit_file <- snakemake@input[["fit"]]

saveRDS(df, file = snakemake@output[["out"]])




