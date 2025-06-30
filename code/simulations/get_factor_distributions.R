library(ieugwasr)
library(dplyr)
library(stringr) 
library(tidyr)
library(ashr)


files <- paste0("../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.", 1:22, ".RDS")
ref_path <- "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR"

# random LD pruning
set.seed(54321)
P <- purrr::map_df(files, function(f){
                       cat(f, "\n")
                       x <- readRDS(f)
                       Pi <- select(x, chrom, snp, pos) %>% rename(rsid = snp)
                       Pi$pval <- runif(n = nrow(Pi), min = 0, max = 1)
                       ld_res <- ld_clump(Pi, plink_bin = genetics.binaRies::get_plink_binary(), bfile = ref_path,
                                          clump_p = 1, clump_r2 = 0.01, clump_kb = 1000)
                       filter(x, snp %in% ld_res$rsid)
                       })
#dim(P)
# [1] 71850    33
nfactor <- 14
ash_dists <- lapply(1:nfactor,  function(i){
                       bh <- P[[paste0("factor", i, ".z")]]
                       ash_res <- ash(betahat = bh, sebetahat = 1, mixcompdist = "normal")$fitted_g
                       mypi <- ash_res$pi[ash_res$pi !=0]
                       mysd <- ash_res$sd[ash_res$pi != 0]
                       mydist <- GWASBrewer::mixnorm_to_scale_fam(sigma = mysd, pi = mypi)
                       })
saveRDS(ash_dists, file = "bc_factor_dists.RDS")

