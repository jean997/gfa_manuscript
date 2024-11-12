library(dplyr)
library(ieugwasr)

X <- readRDS(snakemake@input[["zmat"]])
r2_thresh <- as.numeric(snakemake@wildcards[["r2_thresh"]])
clump_kb <- snakemake@wildcards[["kb"]]
ref_path  <- snakemake@params[["ref_path"]]
out <- snakemake@output[["out"]]
p <- snakemake@wildcards[["p"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])

jitter_sd <- as.numeric(snakemake@wildcards[["jitter_sd"]])
jitter_seed <- as.numeric(snakemake@wildcards[["jitter_seed"]])

set.seed(jitter_seed) # this only matters if we are jittering


if(!p %in% c("pvalue", "rank")){
  stop("Unknown prioritization option.\n")
}


Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()

if(p == "pvalue"){
  zmax <- apply(abs(Z_hat), 1, function(x){max(x, na.rm=T)})
  if(jitter_sd > 0){
    zmax <- zmax + rnorm(n = length(zmax), sd = jitter_sd, mean = 0)
  }
  myp <- 2*pnorm(-abs(zmax))

}else if(p == "rank"){
  Z_rank <- apply(Z_hat,2,function(x){rank(x,na.last = "keep")})
  min_rank <- apply(Z_rank, 1, function(x){min(x, na.rm = T)})
  myp <- min_rank/max(min_rank)

}

dat <- data.frame(rsid = X$snp, pval = myp)



dat_clump <- ld_clump(dat = dat,
                     clump_r2 = r2_thresh,
                     clump_p = pthresh,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_path)

ix <- which(X$snp %in% dat_clump$rsid)
X <- X[ix,]

saveRDS(X, file=out)

