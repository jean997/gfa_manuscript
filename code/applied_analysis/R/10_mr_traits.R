source("renv/activate.R")
library(ieugwasr)
library(MVMR)
library(dplyr)
library(ieugwasr)
library(purrr)
library(stringr)
library(gwasvcf)

source("R/format_ieu_chrom.R")

z_files = unlist(snakemake@input[["Z"]])
r2_thresh <- as.numeric(snakemake@params[["r2_thresh"]])
clump_kb <- snakemake@params[["clump_kb"]]
ref_path  <- snakemake@params[["ref_path"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])
outcome <- snakemake@wildcards[["outcome"]]
outcome_info_file <- snakemake@params[["outcome_info_file"]]
info <- read_csv(outcome_info_file)
R <- readRDS(snakemake@input[["R"]])

out <- snakemake@output[["out"]]

i <- which(info$name == outcome)
outcome_file <- info$raw_data_path[i]

index_traits <- readRDS("index_traits.RDS")
mycols <- c(paste0(index_traits, ".z"), "snp", "chrom")
X <- map_dfr(1:22, function(c){
              f <- z_files[c]
              x <- readRDS(f) %>% dplyr::select(all_of(mycols))

              if(str_ends(outcome_file, ".vcf.gz") | str_ends(outcome_file, ".vcf")){
                out_dat <- format_ieu_chrom(outcome_file, c, 0) %>%
                    dplyr::select(snp, beta_hat, se, p_value) %>% 
                    dplyr::rename(beta_out = beta_hat, se_out = se, p_out = p_value)
              }else{
                out_dat <- format_flat_chrom(outcome_file, c, 0,
                                                     info$snp[i],
                                                     info$pos[i],
                                                     info$chrom[i],
                                                     info$A1[i],
                                                     info$A2[i],
                                                     info$beta_hat[i],
                                                     info$se[i],
                                                     info$p_value[i],
                                                     info$af[i],
                                                     info$sample_size[i],
                                                     as.logical(info$effect_is_or[i]))
                out_dat <- out_dat %>% 
                    dplyr::select(snp, beta_hat, se, p_value) %>% 
                    dplyr::rename(beta_out = beta_hat, se_out = se, p_out = p_value)
              }


              full_dat <- inner_join(out_dat, x, by = "snp")
              z <- dplyr::select(full_dat, ends_with(".z")) %>% as.matrix()
              minp <- 2*pnorm(-abs(apply(abs(z), 1, max)))
              ld_res <- ld_clump(data.frame(rsid = full_dat$snp, pval = minp), 
                                 plink_bin = genetics.binaRies::get_plink_binary(), 
                                 bfile = ref_path, 
                                 clump_p = pthresh, 
                                 clump_r2 = r2_thresh, 
                                 clump_kb = clump_kb)
              return(filter(full_dat, snp %in% ld_res$rsid))
})

nexp <- length(grep(".z$", names(X)))
trait_names <- names(dplyr::select(X, ends_with(".z"))) %>% str_replace(".z$", "")

mvmr_dat <- format_mvmr(BXGs =  (dplyr::select(X, ends_with(".z")) %>% as.matrix()), 
                        BYG = X$beta_out, 
                        seBXGs = matrix(1, nrow = nrow(X), ncol = nexp),
                        seBYG = X$se_out,
                        RSID = X$snp)
trait_names <- names(dplyr::select(X, ends_with(".z"))) %>% str_replace(".z$", "")

o <- match(trait_names, R$names)
Rcor <- Matrix::nearPD(cov2cor(R$R), corr = TRUE)
Rcor <- as.matrix(Rcor$mat)
gencov <- lapply(seq(nrow(X)), function(x){Rcor[o,o]})


str_res <- strength_mvmr(mvmr_dat, gencov = gencov)
ivw_res <- ivw_mvmr(mvmr_dat, gencov = gencov)

ret <- list(mvmr_dat = mvmr_dat, str_res = str_res, ivw_res = ivw_res, trait_names = trait_names, gencov = gencov)

saveRDS(ret, file = out)
