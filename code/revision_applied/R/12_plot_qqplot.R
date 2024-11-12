library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)


# sample size affects genetic covariance and h2 but not intercept or genetic correlation
z_files = unlist(snakemake@input[["Z"]])
qq_out_prefix <- snakemake@output[["qq_out"]] |> str_replace(".1.png$", "")

X <- map_dfr(z_files, function(f){
  readRDS(f) 
})

fct_names <- str_subset(names(X), "\\.p$")
fct_names <- str_replace(fct_names, "\\.p$", "")
nf <- length(fct_names)
# order as shown in publication figures
if(str_detect(qq_out_prefix, "metab_")){
    fct_order <- c(1,  5, 15, 11,  2, 10,  7,  4, 13,  3,  8,  6,  9, 14, 12)
    disp_name <- paste0("Factor ", 1:nf)
}else if(str_detect(qq_out_prefix, "bc_")){
    fct_order <- c(7,  2,  3,  6,  1,  9, 14, 10,  5, 11,  4,  8, 12, 13) #c(7, 2, 1,3,6, 14, 10, 9, 5,11, 4, 8, 12, 13)
    disp_name <- paste0("Factor ", 1:nf)
}


nplot <- ceiling(nf/4)

i <- 1
for(j in 1:nplot){
    png(paste0(qq_out_prefix, ".", j, ".png"), 
        height = 8, width = 8, units = "in", res = 300)
    par(mfrow = c(2, 2))
    nmax <- min(4, length(fct_names)-i -1)
    for(k in 1:nmax){
        qqman::qq(X[[paste0(fct_names[fct_order[i]], ".p")]], main = disp_name[i])
        i <- i + 1
    }
    dev.off()
}

