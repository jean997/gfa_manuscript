source("renv/activate.R")

inp <- unlist(snakemake@input)
inp <- unique(inp)
res <- purrr::map_df(inp, readRDS)
saveRDS(res, file = snakemake@output[["out"]])

