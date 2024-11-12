library(GFA)
library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)


source("R/order_scale_factors.R")

res <- readRDS(snakemake@input[["gfa_fit"]])
names_simple <- str_replace(res$names, "Astle_2016_2786352__", "")
trait_order <- c("basophil-count",
                 "basophil-percentage-of-granulocytes",
                 "basophil-percentage-of-white-cells",
                 "eosinophil-count",
                 "eosinophil-percentage-of-white-cells",
                 "granulocyte-percentage-of-myeloid-white-cells",
                 "sum-eosinophil-basophil-counts",
                 "sum-neutrophil-eosinophil-counts",
                 "neutrophil-count",
                 "neutrophil-percentage-of-granulocytes",
                 "neutrophil-percentage-of-white-cells",
                 "lymphocyte-count",
                 "lymphocyte-percentage-of-white-cells",
                 "monocyte-count",
                 "monocyte-percentage-of-white-cells",
                 "myeloid-white-cell-count",
                 "white-blood-cell-count",
                 "platelet-count",
                 "platelet-distribution-width",
                 "mean-platelet-volume",
                 "plateletcrit" ,
                 "hematocrit",
                 "hemoglobin-concentration",
                 "mean-corpuscular-hemoglobin-concentration",
                 "mean-corpuscular-hemoglobin",
                 "mean-corpuscular-volume" ,
                 "reticulocyte-count",
                 "high-light-scatter-reticulocyte-count" ,
                 "immature-reticulocyte-fraction" ,
                 "red-blood-cell-count",
                 "red-cell-distribution-width"

                                            )

o <- match(trait_order, names_simple)
o <- o[!is.na(o)]
names_simple <- str_replace_all(names_simple, "-", " ")
trait_order <- str_replace_all(trait_order, "-", " ")

abbrv <- function(nms){
nms <- str_replace_all(nms, "counts", "#s")
nms <- str_replace_all(nms, "count", "#")
nms <- str_replace_all(nms, "percentage", "%")
nms <- str_replace_all(nms, "fraction", "frac")
nms <- str_replace_all(nms, "distribution", "dist")
nms <- str_replace_all(nms, "concentration", "conc")
return(nms)
}

names_simple2 <- abbrv(names_simple)

F_hat <- order_and_scale(res$F_hat[o,])
fct_order <- F_hat$order

fct_plot <- F_hat$F_hat %>%
  plot_factors(row_names = names_simple2[o]) +
  scale_fill_gradient2(low = viridis(4)[3], mid = "white", high = viridis(4)[1],
                       name = "Effect") +
  xlab("Factor") + ylab("Trait") +
  theme(axis.text.x  = element_text(size=14, angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size=15),
        strip.text = element_text(size=15))
ggsave(fct_plot, file = snakemake@output[["factor_plot"]], create.dir = TRUE,
       height = 6, width = 9, units = "in", dpi = 300)


if(!res$method == "noR"){
## genetic correlation of factors
gencor_fct <- readRDS(snakemake@input[["fct_gencor"]])

gencor_fct_plot <- gencor_fct$Rg[fct_order, fct_order] %>%
  plot_factors(row_title = "Factor") +
  scale_fill_gradient2(name = "Gen Cor") +
  theme(axis.text.x  =element_text(size=14, angle = 0),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size=15),
        strip.text = element_text(size=15))

ggsave(gencor_fct_plot, file = snakemake@params[["fct_gencor_plot"]],  create.dir = TRUE,
       height = 5, width = 7, units = "in", dpi = 300)



### Overlap of top snps
factor_ld_pruned_files = snakemake@input[["factor_ld_pruned_files"]]

factor_ov <- purrr::map_dfr(factor_ld_pruned_files, function(f){
                                 readRDS(f)$overlap_table}) %>% 
                  group_by(id1, id2) %>% summarize(n = sum(n)) 

factor_tot <- factor_ov %>% filter(id1 == id2) %>% select(id1, id2, n) %>% rename(ntot = n) %>% ungroup()

fov <- factor_ov %>%
  ungroup() %>%
  left_join(select(factor_tot, id1, ntot), by = "id1") %>%
  left_join(select(factor_tot, id2, ntot), by = "id2") %>%
  mutate(id1 = factor(id1, levels = paste0("factor", 1:14)),
         id2 = factor(id2, levels = paste0("factor", 1:14)),
         n = n/sqrt(ntot.x*ntot.y)) %>%
  ggplot(., aes(x = id1, y = id2, fill = n)) +
  geom_tile() +
  scale_fill_viridis(name = "Overlap Percent") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
ggsave(fov, file = snakemake@params[["factor_overlap_plot"]], create.dir = TRUE,
       height = 5, width = 7, units = "in", dpi = 300)

}
