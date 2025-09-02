library(dplyr)
library(GFA)
library(stringr)
library(purrr)

input_files <- unlist(snakemake@input)
jitter <- str_split(input_files, "jitter") %>% map(2) %>% unlist() %>% str_split("_") %>% map(1) %>% unlist()
jitter <- (jitter != "0") + 0

gfaseed <- str_split(input_files, "gfaseed") %>% map(2) %>% unlist() %>% str_split("_") %>% map(1) %>% unlist()
gfaseed <- as.numeric(gfaseed)

res_df <- data.frame(file = input_files, 
                     jitter = jitter, 
                     gfaseed = gfaseed)

ix0 <- with(res_df, which(jitter == 0 & gfaseed == 1))

F0 <- readRDS(res_df$file[ix0])$F_hat
F_hats <- lapply(res_df$file, function(f){ readRDS(f)$F_hat})
res_df$nfactor <- sapply(F_hats, function(Fh){ncol(Fh)})

min_norms <- lapply(F_hats, function(Fh){min_norm(f_true = F0, f_hat = Fh, single_trait_thresh=1)})
res_df$frob_n <- sapply(min_norms, function(x){x$frob_n})
res_df$worst_cor <- sapply(min_norms, function(x){min(x$solution$val[x$solution$val > 0], na.rm=T)})
res_df$n_missed  <- sapply(min_norms, function(x){sum(x$solution$val < 0.8 & !is.na(x$solution$true_ix))})
saveRDS(res_df, file = snakemake@output[["out"]])

