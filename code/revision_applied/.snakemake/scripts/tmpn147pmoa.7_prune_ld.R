
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
    input = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_ldsc.7.RDS', "Z" = '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_ldsc.7.RDS'),
    output = list('../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_ldsc.top_snps.7.RDS', "out" = '../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_ldsc.top_snps.7.RDS'),
    params = list('/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', 0.01, 1000, '5e-8', "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "r2_thresh" = 0.01, "clump_kb" = 1000, "pthresh" = '5e-8'),
    wildcards = list('bc', 'gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_ldsc', '7', "prefix" = 'bc', "analysis" = 'gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_ldsc', "chrom" = '7'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("input" = list("metab" = 'orig_metab_traits.csv', "bc" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt', 'none'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = 100), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = 1, "af_thresh" = 0.01, "sample_size_tol" = 0.1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'topsnps_fctgls',
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
library(dplyr)
library(stringr)
library(tidyr)

dat <- readRDS(snakemake@input[["Z"]])
r2_thresh <- as.numeric(snakemake@params[["r2_thresh"]])
clump_kb <- snakemake@params[["clump_kb"]]
ref_path  <- snakemake@params[["ref_path"]]
out <- snakemake@output[["out"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])

P <- dat %>%
         select(ends_with(".z")) %>%
         mutate_all(~2*pnorm(-abs(.x)))
names(P) <- str_replace(names(P), ".z$", "") 
Pnames <- names(P)
P <- bind_cols(select(dat, "snp"), P)
P <- pivot_longer(data = P, names_to = "id", values_to = "pval", cols = all_of(Pnames)) %>%
     rename(rsid = snp) %>%
     filter(pval < pthresh)

Pu <- P %>% group_by(rsid) %>% summarize(pval = min(pval))

ld_res <- ld_clump(P, plink_bin = genetics.binaRies::get_plink_binary(), bfile = ref_path, 
                   clump_p = pthresh, clump_r2 = r2_thresh, clump_kb = clump_kb)

ld_resu <- ld_clump(Pu, plink_bin = genetics.binaRies::get_plink_binary(), bfile = ref_path, 
                   clump_p = pthresh, clump_r2 = r2_thresh, clump_kb = clump_kb)

P_top <- filter(P, rsid %in% ld_res$rsid)
overlap_table <- purrr::map_dfr(Pnames, function(n1){
    top_snps <- ld_res %>% filter(id == n1) %>% pull(rsid)
    if(length(top_snps) == 0){
        return(data.frame(id1 = n1, id2 = Pnames, n = 0))
    }
    P_top <- filter(P, rsid %in% top_snps)
    nov <-  filter(P_top, pval < pthresh) %>% 
            group_by(id) %>% 
            summarize(n = n())
    names(nov)[1] <- c("id2")
    nov$id1 <- n1
    nov <- nov %>% select(id1, id2, n)
    if(any(!Pnames %in% nov$id2)){
        nms <- Pnames[!Pnames %in% nov$id2]
        nov <- bind_rows(nov, data.frame(id1 = n1, id2 = nms, n = 0))
    }
    return(nov)
})

result <- list(ld_res = ld_res, ld_resu = ld_resu, overlap_table = overlap_table)
saveRDS(result, file = out)
#ggplot(overlap_table, aes(x = id1, y = id2, fill = n)) + geom_tile() + scale_fill_viridis_b()








