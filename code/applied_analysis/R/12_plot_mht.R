source("renv/activate.R")
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(topr)

source("renv/activate.R")
z_files_factors = unlist(snakemake@input[["Z_factors"]])
z_files_traits = unlist(snakemake@input[["Z_traits"]])
res <- readRDS(snakemake@input[["gfa_res"]])
mht_out_prefix <- snakemake@output[["mht_out"]] |> str_replace(".1.png$", "")

X_fct <- map_dfr(z_files_factors, function(f){
  readRDS(f) 
})

X_fct$chrom <- as.numeric(X_fct$chrom)

fct_names <- str_subset(names(X_fct), "\\.p$")
fct_names <- str_replace(fct_names, "\\.p$", "")

nf <- length(fct_names)

# order as shown in publication figures
if(str_detect(mht_out_prefix, "metab_")){
    fct_order <- c(1,  5, 15, 11,  2, 10,  7,  4, 13,  3,  8,  6,  9, 14, 12)
    disp_name <- paste0("CAD/T2D Factor ", 1:nf)
}else if(str_detect(mht_out_prefix, "bc_")){
    fct_order <- c(7,  2,  3,  6,  1,  9, 14, 10,  5, 11,  4,  8, 12, 13) 
    disp_name <- paste0("Blood Cell Factor ", 1:nf)
}

X_trait <- map_dfr(z_files_traits, function(f){
  readRDS(f) 
})
trait_names <- str_subset(names(X_trait), "\\.z$")
trait_names <- str_replace(trait_names, "\\.z$", "")

disp_trait_names <- str_split(trait_names, "__") %>% map(2) %>% unlist() %>% str_replace("-", " ")


#nplot <- ceiling(nf/2)


i <- 1
for(j in 1:nf){
    out_file <- paste0(mht_out_prefix, ".", j, ".png") 
    which_f <- paste0("factor", fct_order[j], ".p")
    lookup <- c(CHROM = "chrom", POS = "pos", "P" = which_f)
    myX <- select(X_fct, all_of(c("chrom", "pos", which_f))) %>%
           rename(all_of(lookup))

    tname <- res$names[which.max(abs(res$F_hat[, fct_order[j]]))]
    which_trait <- paste0(tname, ".z")
    lookup_trait <- c("Z" = which_trait)
    myTrt <- select(X_trait, all_of(c(which_trait))) %>%
             rename(all_of(lookup_trait)) %>%
             mutate(P = 2*pnorm(-abs(Z))) %>%
             select(P)
    myTrt <- cbind(myX[, c("CHROM", "POS")], myTrt)
    trait_index <- match(tname, trait_names)
    disp_tname <- disp_trait_names[trait_index]
    disp_fname <- disp_name[j]
    plt <- manhattan(list(myX, myTrt), annotate = 1e-20, legend_labels = c(disp_fname, disp_tname), ntop =1, title = disp_fname, build = 37)
    ggsave(plt, file= out_file, height = 8, width = 15, units = "in", dpi = 300)
}





