
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
    input = list('/nfs/turbo/sph-jvmorr/gwas_summary_statistics/van-der-Harst_2017_29212778/ebi-a-GCST005195.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Xue_2018_30054458/30054458-GCST006867-EFO_0001360.h.tsv.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-19953.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Liu_2019_30643251/ieu-b-4877.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Liu_2019_30643251/ieu-b-73.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-3957.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-16781.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Evangelou_2018_30224653/ieu-b-38.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Evangelou_2018_30224653/ieu-b-39.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-9405.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-8909.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Richardson_2020_32203549/ieu-b-110.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Richardson_2020_32203549/ieu-b-109.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Richardson_2020_32203549/ieu-b-111.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Said_2022_35459240/GCST90029070_buildGRCh37.tsv.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-15590.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-16446.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-4650.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-8476.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Trubetskoy_2022_35396580/ieu-b-5099.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Demontis_2019_30478444/ieu-a-1183.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Mullins_2021_34002096/ieu-b-5110.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Howard_2019_30718901/ieu-b-102.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-10787.vcf.gz', 'orig_metab_traits.csv', "files" = c('/nfs/turbo/sph-jvmorr/gwas_summary_statistics/van-der-Harst_2017_29212778/ebi-a-GCST005195.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Xue_2018_30054458/30054458-GCST006867-EFO_0001360.h.tsv.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-19953.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Liu_2019_30643251/ieu-b-4877.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Liu_2019_30643251/ieu-b-73.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-3957.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-16781.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Evangelou_2018_30224653/ieu-b-38.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Evangelou_2018_30224653/ieu-b-39.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-9405.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-8909.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Richardson_2020_32203549/ieu-b-110.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Richardson_2020_32203549/ieu-b-109.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Richardson_2020_32203549/ieu-b-111.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Said_2022_35459240/GCST90029070_buildGRCh37.tsv.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-15590.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-16446.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-4650.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-8476.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Trubetskoy_2022_35396580/ieu-b-5099.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Demontis_2019_30478444/ieu-a-1183.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Mullins_2021_34002096/ieu-b-5110.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Howard_2019_30718901/ieu-b-102.vcf.gz', '/nfs/turbo/sph-jvmorr/gwas_summary_statistics/Elsworth_2018_0/ukb-b-10787.vcf.gz'), "gwas_info" = 'orig_metab_traits.csv'),
    output = list('../../data/gfa_intermediate_data/metab_zmat.22.RDS', "out" = '../../data/gfa_intermediate_data/metab_zmat.22.RDS'),
    params = list(0.01, 0.1, "af_thresh" = 0.01, "sample_size_tol" = 0.1),
    wildcards = list('metab', '22', "prefix" = 'metab', "chrom" = '22'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 12367, "mem_mib" = 11795, "disk_mb" = 12367, "disk_mib" = 11795, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = 100), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = 1, "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
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

