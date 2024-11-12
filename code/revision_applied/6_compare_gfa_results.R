library(dplyr)
library(GFA)
library(stringr)
library(purrr)


main_file <- snakemake@input[["main_input"]]
other_files <- snakemake@input[["other_files"]]


seeds <-  str_split(other_files, "gfaseed") %>% map(2) %>% unlist() %>% str_split(fixed("_")) %>% map(1) %>% unlist()
jitter_sd <-  str_split(other_files, "jitter") %>% map(2) %>% unlist() %>% str_split(fixed("_")) %>% map(1) %>% unlist()
jitter_seed <- str_split(other_files, "_") %>% map(11) %>% str_replace(".R$", "")
cond_num <-  str_split(other_files, "cn") %>% map(2) %>% unlist() %>% str_split(fixed(".")) %>% map(1) %>% unlist()


res <- data.frame(gfa_seed = seeds, jitter_sd = jitter_sd, jitter_seed = jitter_seed, cond_num = cond_num)

main_res <- readRDS(main_file)

solutions  <- lapply(other_files, function(o){
                         other_res <- readRDS(o)
                         sol <- min_norm(f_hat = other_res$F_hat, f_true = main_res$F_hat)
                         sol})
res$frob_n <- map(solutions, "frob_n")




