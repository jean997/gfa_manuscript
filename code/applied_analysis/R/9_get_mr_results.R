library(stringr)
library(dplyr)
library(viridis)

mr_files_factors <- snakemake@input[["mr_input_factors"]]
mr_files_traits <- snakemake@input[["mr_input_factors"]]
res <- readRDS(snakemake@input[["gfa_fit"]]) # needed to scale factor effects to be comparable to trait effects

fct_order <- c(1, 2, 4, 7, 9, 8, 10, 14, 3, 5, 6, 11, 12, 13)
# names of index traits
ix_names2 <- c("eosinophil #","neutrophil #","lymphocyte #","monocyte #",
 "basophil #","platelet #", "mean platelet volume","platelet dist width",
 "reticulocyte #","mean corpuscular hemoglobin", "hematocrit",
 "mean corpuscular hemoglobin conc", "immature reticulocyte frac","red cell dist width")     

ix_names1 <- c("eosinophil-count","neutrophil-count","lymphocyte-count","monocyte-count",
 "_basophil-count","platelet-count", "mean-platelet-volume","platelet-distribution-width",
 "_reticulocyte-count","mean-corpuscular-hemoglobin$", "hematocrit",
 "mean-corpuscular-hemoglobin-concentration", "immature-reticulocyte-fraction","red-cell-distribution-width")     

res_index_trait_location <- sapply(ix_names1, function(n){grep(n, res$names)}) %>% as.numeric()

fct_dict <- data.frame(factor = paste0("factor", fct_order),
                       index_trait = ix_names2,
                       plot_order = 1:14)
fct_dict$fct_scale_factor <- sapply(1:14, function(j){
    f <- fct_order[j]
    v <- res_index_trait_location[j]
    sf <- as.numeric(res$F_hat[v, f])
    return(sf)
  })



outcome_traits <- c("type-1-diabetes", "asthma", "celiac-disease", "inflammatory-bowel-disease", "rheumatoid-arthritis", 
                    "multiple-sclerosis", "CAD", "type-2-diabetes")

mr_files_factors <- mr_files_factors[ as.numeric(sapply(outcome_traits, function(o){grep(o, mr_files_factors)}))]
mr_files_traits <- mr_files_traits[ as.numeric(sapply(outcome_traits, function(o){grep(o, mr_files_traits)}))]

all_mr_results <- purrr::map_dfr(seq_along(outcome_traits), function(i){
  res_fct <- readRDS(mr_files_factors[i])
  res_trait <- readRDS(mr_files_traits[i])

  mr_results <- fct_dict %>%
    mutate(outcome = outcome_traits[i])

  mr_results$fct_beta <- res_fct$ivw_res[str_replace(mr_results$factor, "factor", "exposure"), "Estimate"]*mr_results$fct_scale_factor
  mr_results$fct_se <- res_fct$ivw_res[str_replace(mr_results$factor, "factor", "exposure"), "Std. Error"]*mr_results$fct_scale_factor
  mr_results$fct_p <- res_fct$ivw_res[str_replace(mr_results$factor, "factor", "exposure"), "Pr(>|t|)"]

  mr_results$trait_beta <- res_trait$ivw_res[, "Estimate"]
  mr_results$trait_se <- res_trait$ivw_res[, "Std. Error"]
  mr_results$trait_p <- res_trait$ivw_res[, "Pr(>|t|)"]
  mr_results$trait_F <- unlist(res_trait$str_res)
  mr_results$fct_F <- unlist(res_fct$str_res)
  return(mr_results)
})


pub_table <- select(all_mr_results, plot_order, index_trait, outcome, fct_beta, fct_se, fct_p, trait_beta, trait_se, trait_p, trait_F, fct_F) %>%
  rename(Factor=plot_order,
         `Index Trait` = index_trait,
         Outcome = outcome,
         `Factor Beta` = fct_beta,
         `Factor SE` = fct_se,
         `Factor P` = fct_p,
         `Trait Beta` = trait_beta,
         `Trait SE` = trait_se,
         `Trait P` = trait_p,
         `Trait F-stat` = trait_F,
         `Factor F-stat`  = fct_F)
readr::write_csv(pub_table, file = snakemake@output[["out"]])
