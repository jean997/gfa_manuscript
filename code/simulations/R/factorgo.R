source("renv/activate.R")
library(GWASBrewer)
library(readr)


dat <- readRDS(snakemake@input[["inp"]])
pthresh <- as.numeric(snakemake@wildcards[["pthresh"]])
nmax <- as.numeric(snakemake@params[["nmax"]])
minsig <- as.numeric(snakemake@wildcards[["minsig"]])




###

r2_thresh <- 0.01
n <- ncol(dat$beta_hat)
o <- c(n + 1, 1:n)

fgo_prefix <- snakemake@params[["prefix"]]
dat_dir <- snakemake@params[["temp_dir"]]

rand_string <- function(n) {
  paste0(sample(c(letters, LETTERS, 0:9), n, replace = TRUE), collapse = "")
}
dat_prefix <- paste0(dat_dir, rand_string(10))

zhat_file <- paste0(dat_prefix, ".zhat.tsv.gz")
N_file <- paste0(dat_prefix, ".N.tsv")


# for factorGo, use only SNPs significant for at least {minsig} traits
# select these before LD pruning
pmin <- apply(dat$pval, 1, min);
nsig2 <- rowSums(dat$pval < pthresh)
pmin[nsig2 < minsig] <- 1
ix <- sim_ld_prune(dat, pvalue = pmin,
                    R_LD = ld_mat_list,
                    r2_thresh = r2_thresh,
                    pval_thresh = pthresh);

Z_hat <- with(dat, beta_hat[ix,]/s_estimate[ix,])
Z_hat <- data.frame(Z_hat)
Z_hat$rsid <- ix
Z_hat <- Z_hat[, o]
readr::write_tsv(Z_hat, file = zhat_file)


N <- data.frame("N" = dat$N)
readr::write_tsv(N, file = N_file)

nmax <- min(nmax, floor(nrow(Z_hat)/3))


#get list of nF
s <- svd(Z_hat);
nF <- try(nFactors::nScree(x = s$d),
          silent = TRUE)
# this deals with one-off error in nScree that occurs when nparallel = n
if(inherits(nF, "try-error")){
    nF <- nFactors::nScree(x = c(s$d, 0))
}


nF <- as.list(nF$Components)
i <- which(names(nF) == "nparallel")
nF <- nF[-i]
nF$six <- 6
nF$twelve <- 12

fgo_init_prefix <- paste0(fgo_prefix, ".init")
fgo_cmd1 <- paste0("factorgo ", zhat_file, " ", N_file, " -k  ", nmax, " --scale -o  ", fgo_init_prefix)
system(fgo_cmd1)
fct <- read_table(paste0(fgo_init_prefix, ".factor.tsv.gz"), col_names = FALSE)
nguess <- sum(fct$X3 > 0.001)
system(paste0("rm ", fgo_init_prefix, "*"))

nF$guess <- nguess

fits <- lapply(nF, function(n){
    if(n == 0){
        fit <- list("F_hat" = matrix(NA, nrow = ncol(Z_hat), ncol = 0), 
                    "L_hat" = NULL,
                     "ix" = ix) 
        return(fit)
    }
    fgo_cmd2 <- paste0("factorgo ", zhat_file, " ", N_file, " -k  ", min(n, nmax), " --scale -o  ", fgo_prefix)
    system(fgo_cmd2)
    F_hat <- read_table(paste0(fgo_prefix, ".Zm.tsv.gz"), col_names = F) |> as.matrix()
    L_hat <- read_table(paste0(fgo_prefix, ".Wm.tsv.gz"), col_names = F) |> as.matrix()
    F_hat_var <- read_table(paste0(fgo_prefix, ".Zvar.tsv.gz"), col_names = F) |> as.matrix()
    L_hat_var <- read_table(paste0(fgo_prefix, ".Wvar.tsv.gz"), col_names = F) |> as.matrix()
    fct <- read_table(paste0(fgo_prefix, ".factor.tsv.gz"), col_names = F) 

    fit <- list("F_hat" = F_hat, "L_hat" = L_hat,
            "F_hat_var" = F_hat_var, "L_hat_var" = L_hat_var,
            "fct" = fct,
             "ix" = ix)
    fit
                    })


system(paste0("rm ", zhat_file))
system(paste0("rm ", N_file))

saveRDS(fits, file = snakemake@output[["out"]])



