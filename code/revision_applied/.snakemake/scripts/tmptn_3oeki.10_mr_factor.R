
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.22.RDS', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Xue_2018_30054458/30054458-GCST006867-EFO_0001360.h.tsv.gz', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.Rgcor.RDS', "Z" = c('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.1.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.2.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.3.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.4.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.5.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.6.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.7.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.8.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.9.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.10.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.11.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.12.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.13.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.14.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.15.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.16.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.17.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.18.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.19.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.20.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.21.RDS', '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.22.RDS'), "outcome_file" = '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Xue_2018_30054458/30054458-GCST006867-EFO_0001360.h.tsv.gz', "R" = '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.Rgcor.RDS'),
    output = list('../../results/applied_analysis/mr_results/bc_mr_fct.Xue_2018_30054458__type-2-diabetes.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.RDS', "out" = '../../results/applied_analysis/mr_results/bc_mr_fct.Xue_2018_30054458__type-2-diabetes.gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc.RDS'),
    params = list('/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', 5e-08, 0.01, 1000, 'orig_autoimmune.csv', "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "pthresh" = 5e-08, "r2_thresh" = 0.01, "clump_kb" = 1000, "outcome_info_file" = 'orig_autoimmune.csv'),
    wildcards = list('bc', 'Xue_2018_30054458__type-2-diabetes', 'gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc', "prefix" = 'bc', "outcome" = 'Xue_2018_30054458__type-2-diabetes', "analysis" = 'gfaseed1_cn100.ldpruned_r20.01_kb1000_pvalue_jitter0_0.R_ldsc'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 3497, "mem_mib" = 3335, "disk_mb" = 3497, "disk_mib" = 3335, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "jitter_sd" = c(0, 0.5), "jitter_seed" = c(1, 2, 3)), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = c(100)), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1, 2, 3), "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'mr_factor',
    bench_iteration = as.numeric(NA),
    scriptdir = '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_analysis/gfa_manuscript/code/revision_applied/R',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
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


X <- map_dfr(1:22, function(c){
              f <- z_files[c]
              x <- readRDS(f)
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
              return(dplyr::filter(full_dat, snp %in% ld_res$rsid))
})

nexp <- length(grep(".z$", names(X)))
mvmr_dat <- format_mvmr(BXGs =  (dplyr::select(X, ends_with(".z")) %>% as.matrix()), 
                        BYG = X$beta_out, 
                        seBXGs = matrix(1, nrow = nrow(X), ncol = nexp),
                        seBYG = X$se_out,
                        RSID = X$snp)

Rcor <- Matrix::nearPD(cov2cor(R$R), corr = TRUE)
Rcor <- as.matrix(Rcor$mat)
gencov <- lapply(seq(nrow(X)), function(x){Rcor})

str_res <- strength_mvmr(mvmr_dat, gencov = gencov)
ivw_res <- ivw_mvmr(mvmr_dat, gencov = gencov)

ret <- list(mvmr_dat = mvmr_dat, str_res = str_res, ivw_res = ivw_res, gencov = gencov)

saveRDS(ret, file = out)
