library(GFA)
library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)

### Overlap of top snps

abbrv <- function(nms){
nms <- str_replace_all(nms, "counts", "#s")
nms <- str_replace_all(nms, "count", "#")
nms <- str_replace_all(nms, "percentage", "%")
nms <- str_replace_all(nms, "fraction", "frac")
nms <- str_replace_all(nms, "distribution", "dist")
nms <- str_replace_all(nms, "concentration", "conc")
return(nms)
}

ix_names2 <- c("eosinophil #","neutrophil #","lymphocyte #","monocyte #",
 "basophil #","platelet #", "mean platelet volume","platelet dist width",
 "reticulocyte #","mean corpuscular hemoglobin", "hematocrit",
 "mean corpuscular hemoglobin conc", "immature reticulocyte frac","red cell dist width")     


R_ldsc <- readRDS(snakemake@input[["trait_gencor"]])
R_ldsc$names <- str_replace(R_ldsc$names, "Astle_2016_2786352__", "") %>% str_replace_all("-", " ") %>% abbrv()

to2 <- match(ix_names2, R_ldsc$names)
gencor_idx_trait2 <- cov2cor(R_ldsc$Rg[to2, to2])

gencor_plot_ix2 <- gencor_idx_trait2 %>%
  plot_factors(col_title = "Index Trait", row_title = "Index Trait") +
  scale_fill_gradient2(name = "Gen Cor") +
  theme(axis.text.x  = element_text(size=14, angle = 0),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size=15),
        strip.text = element_text(size=15))

ggsave(gencor_plot_ix2, file = snakemake@output[["trait_gencor_plot"]], create.dir = TRUE,
       height = 5, width = 7, units = "in", dpi = 300)





trait_ld_pruned_files = snakemake@input[["trait_ld_pruned_files"]]

trait_ov <- purrr::map_dfr(trait_ld_pruned_files, function(f){
                                 readRDS(f)$overlap_table}) %>% 
                  group_by(id1, id2) %>% summarize(n = sum(n)) 

trait_ov$id1 <- str_replace(trait_ov$id1, "Astle_2016_2786352__", "") %>% str_replace_all("-", " ") %>% abbrv()
trait_ov$id2 <- str_replace(trait_ov$id2, "Astle_2016_2786352__", "") %>% str_replace_all("-", " ") %>% abbrv()

trait_tot <- trait_ov %>% filter(id1 == id2) %>% select(id1, id2, n) %>% rename(ntot = n) %>% ungroup()

plot_trait_overlap <- trait_ov %>%
  ungroup() %>%
  left_join(select(trait_tot, id1, ntot), by = "id1") %>%
  left_join(select(trait_tot, id2, ntot), by = "id2") %>%
  filter(id1 %in% ix_names2 & id2 %in% ix_names2) %>%
         mutate(id1 = factor(id1, levels = ix_names2),
                id2 = factor(id2, levels = ix_names2),
                n = n/sqrt(ntot.x*ntot.y)) %>%
          ggplot(., aes(x = id1, y = id2, fill = n)) +
          geom_tile() +
          scale_fill_viridis(name = "Overlap Percent") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title = element_blank())

ggsave(plot_trait_overlap, file = snakemake@output[["trait_overlap_plot"]], create.dir=TRUE, height = 6.2, width = 8.6, units = "in", dpi = 300)

