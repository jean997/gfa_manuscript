
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
    input = list('/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004603.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/from_midway/bloodcell_2016/pdw_narrow_form.tsv.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004599.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004604.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004628.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004630.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004622.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004618.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004606.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004627.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004625.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004629.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004600.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004601.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004602.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004605.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004607.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004608.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004609.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004610.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004611.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004612.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004613.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004614.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004615.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004617.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004619.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004620.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004621.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004623.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004624.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004626.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004631.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004632.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004633.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004634.vcf.gz', 'orig_bloodcell.csv', "files" = c('/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004603.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/from_midway/bloodcell_2016/pdw_narrow_form.tsv.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004599.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004604.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004628.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004630.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004622.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004618.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004606.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004627.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004625.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004629.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004600.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004601.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004602.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004605.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004607.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004608.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004609.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004610.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004611.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004612.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004613.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004614.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004615.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004617.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004619.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004620.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004621.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004623.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004624.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004626.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004631.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004632.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004633.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Astle_2016_27863252/ebi-a-GCST004634.vcf.gz'), "gwas_info" = 'orig_bloodcell.csv'),
    output = list('../../data/gfa_intermediate_data/bc1_zmat.4.RDS', "out" = '../../data/gfa_intermediate_data/bc1_zmat.4.RDS'),
    params = list(0.01, 0.1, "af_thresh" = 0.01, "sample_size_tol" = 0.1),
    wildcards = list('bc1', '4', "prefix" = 'bc1', "chrom" = '4'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 62772, "mem_mib" = 59865, "disk_mb" = 62772, "disk_mib" = 59865, "tmpdir" = '/tmp'),
    config = list("input" = list("bc1" = 'orig_bloodcell.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = '1e3'), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = 1, "af_thresh" = 0.01, "sample_size_tol" = 0.1, "cor_clust" = 1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/gfa_results/')),
    rule = 'snp_table_chrom',
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
library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(rlang)
library(readr)
library(purrr)
library(stringr)

source("R/format_ieu_chrom.R")

c <- as.numeric(snakemake@wildcards[["chrom"]])
gwas_info_file <- snakemake@input[["gwas_info"]]

af_thresh <- as.numeric(snakemake@params[["af_thresh"]])
sample_size_tol <- as.numeric(snakemake@params[["sample_size_tol"]])
out <- snakemake@output[["out"]]


info <- read_csv(gwas_info_file)
if(!"af" %in% names(info)){
    info$af <- NA
}

fulldat <- map(seq(nrow(info)),   function(i){
                        f <- info$raw_data_path[i]
                        if(str_ends(f, "vcf.gz") | str_ends(f, "vcf.bgz")){
                            dat <- format_ieu_chrom(f, c, af_thresh)
                        }else{
                            dat <- format_flat_chrom(f, c, af_thresh,
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
                        }
                        if(all(is.na(dat$sample_size))){
                           dat$sample_size <- info$pub_sample_size[i]
                        }

                        if(is.finite(sample_size_tol)){
                           m <- median(dat$sample_size)
                           dat <- filter(dat, sample_size > (1-sample_size_tol)*m & sample_size < (1 + sample_size_tol)*m)
                        }
                        n <- info$name[i]
                        se_name <- as_name(paste0(n, ".se"))
                        z_name <- as_name(paste0(n, ".z"))
                        ss_name <- as_name(paste0(n, ".ss"))
                        af_name <- as_name(paste0(n, ".af"))

                        dat$sample_size[is.na(dat$sample_size)] <- as.numeric(info$pub_sample_size)
                        dat <-dat %>%  dplyr::mutate(Z = beta_hat/se) %>%
                               dplyr::rename(REF = A2, ALT = A1) %>%
                               dplyr::select(chrom, snp, REF, ALT,
                                              !!z_name := Z,
                                              !!se_name := se,
                                              !!ss_name := sample_size,
                                              !!af_name := allele_freq)
                 }) %>%
       purrr::reduce(full_join, by = c("chrom", "snp", "REF", "ALT"))

dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
if(length(dup_snps) > 0){
    fulldat <- filter(fulldat, !snp %in% dup_snps)
}



# Table of how traits are missing each SNP for LD clumping
miss <- fulldat %>%
        dplyr::select(ends_with(".z")) %>%
        is.na(.) %>%
        rowSums(.)

#nmiss <- data.frame(snp = fulldat$snp, miss = miss)
#ix <- which(miss <= nmiss_thresh)
ix <- which(miss == 0)

saveRDS(fulldat[ix,], file=out)

#saveRDS(nmiss[ix,], file=nm_out)

