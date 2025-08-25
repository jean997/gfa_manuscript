

dat <- readRDS(snakemake@input[["dat"]])

nsig <- dat$nsig_per_trait
af <- dat$snp_info$AF
h2_explained <- (dat$beta_joint^2)*af*(1-af)*2

quant_df <- purrr::map_dfc(1:ncol(h2_explained), function(i){
                       x <- h2_explained[,i]
                       u <- x[x != 0]
                       qs <- quantile(u, probs = seq(0.1, 1, by = 0.05))
                       d = data.frame(x = as.numeric(qs))
                       names(d) <- paste0("x", i)
                       return(d)})
quant_df$quantile <- seq(0.1, 1, by = 0.05)


res <- list(nsig = nsig, quant_df = quant_df)

saveRDS(res, file = snakemake@output[["out"]])



