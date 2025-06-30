library(dplyr)
library(GFA)

eval <- function(fit, dat, n_subset = NULL){
    disc_thresh = c(0.9, 0.95, 0.98);
    extra_thresh = c(0.8, 0.7);

    if(!is.null(n_subset)){
        if(n_subset == 0){
            if(!is.null(fit$F_hat)){
                fit$F_hat <- fit$F_hat[, c(), drop = F]
            }
            if(!is.null(fit$L_hat)){
                fit$L_hat <- fit$L_hat[, c(), drop = F]
            }
            fit$B_hat <- NULL
        }else{
            if(!is.null(fit$F_hat)){
                fit$F_hat <- fit$F_hat[, 1:n_subset, drop = F]
            }
            if(!is.null(fit$L_hat)){
                fit$L_hat <- fit$L_hat[, 1:n_subset, drop = F]
            }
            if(!is.null(fit$F_hat) & !is.null(fit$L_hat)){
                fit$B_hat <- fit$L_hat %*% t(fit$F_hat)
            }
        }
    }

    if(is.null(fit$F_hat)){
        disc <- NULL;
        n_disc <- rep(NA, length(disc_thresh));
        n_extra <- rep(NA, length(extra_thresh));
        frob_n <- NA 
        #frob_n_l <- NA
        #err <- NA;
        n_hat <- 0
    }else{
        if(!"matrix" %in% class(fit$F_hat)){
            fit$F_hat <- matrix(fit$F_hat, ncol=1)
            if(!is.null(fit$L_hat)) fit$L_hat <- matrix(fit$L_hat, ncol=1);
        }
        if(is.null(fit$L_hat)){
            err <- NA;
            sol <- min_norm(f_true = dat$F_mat, f_hat = fit$F_hat);
            multi_ix <- sol$solution %>% filter(!is.na(est_ix) & !is.na(val)) %>% pull(est_ix)
            n_hat <- length(multi_ix)
            #frob_n_l <- NA;
        }else{
            sol <- min_norm(f_true = dat$F_mat, f_hat = fit$F_hat,
                                l_true = dat$L_mat_marg[fit$ix,],
                                l_hat = fit$L_hat);
            if(is.null(sol$solution)){
                multi_ix <- NULL
                n_hat <- 0
            }else{
                multi_ix <- sol$solution %>% filter(!is.na(est_ix) & !is.na(val)) %>% pull(est_ix)
                n_hat <- length(multi_ix)
            }
            #frob_n_l <- sol$frob_n_l

            #ix <- fit$ix;
            #Bmulti <- dat$beta_marg[ix,] - dat$theta_marg[ix,];
            #Ball <- dat$beta_marg[ix,]
            #if(is.null(fit$scale)){
            #    B_hat <- (fit$L_hat %*% t(fit$F_hat))*dat$s_estimate[ix,];
            #}else{
            #    B_hat <- (fit$L_hat %*% t(fit$F_hat*fit$scale))*dat$s_estimate[ix,];
            #}
            #err_all <- sqrt(sum((B_hat-Ball)^2)/sum(Ball^2));
    
            #if(length(multi_ix) > 0){
            #    if(is.null(fit$scale)){
            #        B_hatm <- (fit$L_hat[,multi_ix] %*% t(fit$F_hat[,multi_ix]))*dat$s_estimate[ix,];
            #    }else{
            #        B_hatm <- (fit$L_hat[,multi_ix] %*% t(fit$F_hat[,multi_ix]*fit$scale))*dat$s_estimate[ix,];
            #    }
            #    err_multi <- sqrt(sum((B_hatm-Bmulti)^2)/sum(Bmulti^2));
            #}

        }
        n_disc <- sapply(disc_thresh, function(x){sum(sol$solution$val > x, na.rm=T)});
        n_extra <- sapply(extra_thresh, function(x){sum(sol$solution$val < x & !is.na(sol$solution$est_ix), na.rm=T)});
        frob_n <- sol$frob_n
    }
    res <- data.frame(frob_n = frob_n, 
                      n_extra_0.8 = n_extra[1],
                      n_extra_0.7 = n_extra[2],
                      n_disc_0.9  = n_disc[1],                 
                      n_disc_0.95  = n_disc[2],
                      n_disc_0.98   = n_disc[3],
                      n_multi_factors = n_hat)
    return(res)
}

