library(GFA)


res <- readRDS(snakemake@input[["gfa_fit"]])
nsamp <- as.nmumeric(snakemake@params[["nsamp"]])
level <- as.numeric(snakemake@params[["level"]])
pve_int <- GFA:::pve_credints(res, nsamp = nsamp, level = level, type = "pve")
saveRDS(pve_int, file = snakemake@output[["out"]]



