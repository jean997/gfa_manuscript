
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
    input = list('../../data/gfa_intermediate_data/bc97v2_zmat.15.RDS', '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR.bed', "zmat" = '../../data/gfa_intermediate_data/bc97v2_zmat.15.RDS', "bfile" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR.bed'),
    output = list('../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.15.RDS', "out" = '../../data/gfa_intermediate_data/bc97v2_zmat.ldpruned_r20.01_kb1000_pvalue.15.RDS'),
    params = list('/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', 1, 0, "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR', "pthresh" = 1, "is_mvmr" = 0),
    wildcards = list('bc97v2', '0.01', '1000', 'pvalue', '15', "prefix" = 'bc97v2', "r2_thresh" = '0.01', "kb" = '1000', "p" = 'pvalue', "chrom" = '15'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 2270, "mem_mib" = 2165, "disk_mb" = 2270, "disk_mib" = 2165, "tmpdir" = '/tmp'),
    config = list("input" = list("bc95v2" = 'orig_bloodcell_ldsc95.csv', "bc97v2" = 'orig_bloodcell_ldsc97.csv'), "analysis" = list("ldprune" = list("r2_thresh" = 0.01, "clump_kb" = 1000, "ld_prioritization" = c('pvalue'), "ref_path" = '/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR'), "R" = list("type" = c('ldsc', 'pt'), "pthresh" = c(0.05), "l2_dir" = '/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/', "cond_num" = 100), "gfa" = list("gfa_params" = 'default', "max_snps" = 'Inf', "maxrep" = 5, "gfa_seed" = c(1, 2, 3), "af_thresh" = 0.01, "sample_size_tol" = 0.1, "cor_clust" = 1)), "out" = list("data_dir" = '../../data/gfa_intermediate_data/', "output_dir" = '../../results/applied_analysis/gfa_results/')),
    rule = 'ld_prune_plink',
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
library(dplyr)
library(ieugwasr)

X <- readRDS(snakemake@input[["zmat"]])
r2_thresh <- as.numeric(snakemake@wildcards[["r2_thresh"]])
clump_kb <- snakemake@wildcards[["kb"]]
ref_path  <- snakemake@params[["ref_path"]]
out <- snakemake@output[["out"]]
p <- snakemake@wildcards[["p"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])
is_mvmr <- as.numeric(snakemake@params[["is_mvmr"]])

if(!p %in% c("pvalue", "rank")){
  stop("Unknown prioritization option.\n")
}


Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()

if(is_mvmr == 1){
  Z_hat <- Z_hat[,-1,drop = FALSE]
}

if(p == "pvalue"){
  zmax <- apply(abs(Z_hat), 1, function(x){max(x, na.rm=T)})
  myp <- 2*pnorm(-abs(zmax))

}else if(p == "rank"){
  Z_rank <- apply(Z_hat,2,function(x){rank(x,na.last = "keep")})
  min_rank <- apply(Z_rank, 1, function(x){min(x, na.rm = T)})
  myp <- min_rank/max(min_rank)

}

dat <- data.frame(rsid = X$snp, pval = myp)



dat_clump <- ld_clump(dat = dat,
                     clump_r2 = r2_thresh,
                     clump_p = pthresh,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_path)

ix <- which(X$snp %in% dat_clump$rsid)
X <- X[ix,]

saveRDS(X, file=out)

