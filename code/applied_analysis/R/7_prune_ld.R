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








