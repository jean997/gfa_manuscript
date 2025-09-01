library(dplyr)
library(GFA)
library(stringr)
library(purrr)



#### Comparison of results with different seeds and variant sets.
bc_files <- list.files("../../results/applied_analysis/gfa_results/", "bc_gfa", full.names=T) 
bc_files <- stringr::str_subset(bc_files, "ldsc.final.RDS$")
bc_files
# [1] "results/applied_analysis/gfa_results//bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS"  
# [2] "results/applied_analysis/gfa_results//bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_ldsc.final.RDS"
# [3] "results/applied_analysis/gfa_results//bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_2.R_ldsc.final.RDS"
# [4] "results/applied_analysis/gfa_results//bc_gfa_gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_3.R_ldsc.final.RDS"
# [5] "results/applied_analysis/gfa_results//bc_gfa_gfaseed2_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS"  
# [6] "results/applied_analysis/gfa_results//bc_gfa_gfaseed2_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_ldsc.final.RDS"
# [7] "results/applied_analysis/gfa_results//bc_gfa_gfaseed2_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_2.R_ldsc.final.RDS"
# [8] "results/applied_analysis/gfa_results//bc_gfa_gfaseed2_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_3.R_ldsc.final.RDS"
# [9] "results/applied_analysis/gfa_results//bc_gfa_gfaseed3_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.final.RDS"  
# [10] "results/applied_analysis/gfa_results//bc_gfa_gfaseed3_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_1.R_ldsc.final.RDS"
# [11] "results/applied_analysis/gfa_results//bc_gfa_gfaseed3_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_2.R_ldsc.final.RDS"
# [12] "results/applied_analysis/gfa_results//bc_gfa_gfaseed3_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0.5_3.R_ldsc.final.RDS"

F0 <- readRDS(bc_files[1])$F_hat
F_hats <- lapply(bc_files[-1], funciton(f){ readRDS(f)$F_hat})
res_df <- data.frame(gfaseed = rep(1:3, each = 4), jitter = rep(c(0, 1), c(1, 3)))
res_df$nfactor <- c(ncol(F0), sapply(F_hats, function(Fh){ncol(Fh)}))
min_norms <- lapply(F_hats, function(Fh){min_norm(f_true = F0, f_hat = Fh)})
res_df$frob_n <- c(0, sapply(min_norms, function(x){x$frob_n}))
res_df$worst_cor <- c(1, sapply(min_norms, function(x){min(x$solution$val, na.rm=T)}))
saveRDS(res_df, file = paste0(dir, "/gfa_robustness_comparison.RDS"))
#res_df
# gfaseed jitter nfactor       frob_n worst_cor
# 1        1      0      14 0.000000e+00 1.0000000
# 2        1      1      14 1.648637e-01 0.9935537
# 3        1      1      14 1.759755e-01 0.9953691
# 4        1      1      14 1.952024e-01 0.9909897
# 5        2      0      14 3.862841e-16 1.0000000
# 6        2      1      14 1.648637e-01 0.9935537
# 7        2      1      14 1.759755e-01 0.9953691
# 8        2      1      14 1.952024e-01 0.9909897
# 9        3      0      14 3.862841e-16 1.0000000
# 10       3      1      14 1.648637e-01 0.9935537
# 11       3      1      14 1.759755e-01 0.9953691
# 12       3      1      14 1.952024e-01 0.9909897



