source("renv/activate.R")
library(dplyr)
library(ggplot2)
library(tidyr)



res <- readRDS(snakemake@input[["sim_results"]])
dir <- snakemake@params[["plot_dir"]]


res_df_long <- res %>%
  mutate(method = paste0(method, ":", subset),
         scenario = case_when(scenario == "bc" ~ "C", TRUE ~ scenario)) %>%
  select(-data_file, -fit_file, -subset) %>%
  tidyr::pivot_longer(cols = c(paste0("n_disc_", c(0.9, 0.95, 0.98)), paste0("n_extra_", c(0.8, 0.7)), "frob_n", "n_multi_factors"),
                      values_to = "eval", names_to = "eval_type")


## naf for main results
# Figures S7 (j = 0), Main Fig 2 (j = 3), S8 (j = 12)
for(j in c(0, 3, 12)){
  if(j == 3){
    scenarios <-  c("A", "B1", "B2", "B3", "C")
  }else{
    scenarios <- c("A", "B1", "B2", "B3")
  }
  p_res1 <- res_df_long %>%
            filter(method %in% c("gfa-none_pval_1_cn100:all", "gfa-ldsc_pval_1_cn100:all", "gfa-oracle_pval_1_cn100:all",
                                  "ml-all_cn100:naf", "guide-0_5e-8:naf",  "factorgo-2_5e-8:naf", #"factorgo-1_5e-8:naf",
                                  "spc-5e-8:naf", "svd-0_5e-8:naf", "svd-4_5e-8:naf"),
                    eval_type %in% c("n_disc_0.9", "n_extra_0.7", "frob_n"),
                    scenario %in% scenarios,
                    ndense == j | scenario == "C") %>% #, "eval.rrmse")) %>%
            mutate(eval_type = forcats::fct_recode(eval_type,  "Discovered Factors" = "n_disc_0.9",
                                                   "Extraneous Factors" = "n_extra_0.7",
                                                   "Mismatch Score" = "frob_n"),
                   eval_type = factor(eval_type, levels = c("Discovered Factors", "Extraneous Factors", "Mismatch Score")), # "RRMSE")),
                   Method = forcats::fct_recode(method,
                                                "GFA-Rldsc" = "gfa-ldsc_pval_1_cn100:all",
                                                "GFA-Roracle" = "gfa-oracle_pval_1_cn100:all",
                                                "EBMF" = "gfa-none_pval_1_cn100:all",
                                                "MLFA" =  "ml-all_cn100:naf",
                                                "SPC" = "spc-5e-8:naf",
                                                "SVD" = "svd-0_5e-8:naf",
                                                "SVD-HT" = "svd-4_5e-8:naf",
                                                "Guide" = "guide-0_5e-8:naf",
                                                #"FactorGO-1" = "factorgo-1_5e-8:naf",
                                                "FactorGo" = "factorgo-2_5e-8:naf"),
                   Method = factor(Method, levels = c("GFA-Rldsc",  "GFA-Roracle",
                                                      "EBMF", "MLFA", "FactorGo", #,  "FactorGO-1",
                                                      "SPC", "SVD", "SVD-HT", "Guide"))) %>%
            ggplot(aes(x = scenario, y = eval,   color = Method), position = "dodge") +
            geom_boxplot(outlier.size = 1) +
            #geom_violin() +
            xlab("Scenario") +
            ylab("          Score           Num. Extra       Num. Discovered")+
            scale_color_manual( values = as.vector(palette.colors(palette = "Tableau 10", n = 10))) +
            theme_bw() +
            theme(axis.text.y = element_text(size = 10),
                  axis.text.x = element_text(size = 20),
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.title = element_blank(),
                  strip.text = element_text(size = 15)) +
            facet_wrap(~eval_type, ncol = 1, scale = "free_y",
            )

  ggsave(p_res1, file = paste0(dir, "/main_sims_ndense", j, ".png"), width = 12, height = 5, units = "in", dpi = 300)
}




## Null Scenarios
## Figure S9
p_null <- res_df_long %>%
  filter(method %in% c("gfa-none_pval_1_cn100:all", "gfa-ldsc_pval_1_cn100:all", "gfa-oracle_pval_1_cn100:all",
                       "ml-all_cn100:naf", "guide-0_5e-8:naf",
                       "spc-5e-8:naf", "svd-0_5e-8:naf", "svd-4_5e-8:naf"),
         eval_type %in% c("n_extra_0.7", "frob_n"),
         scenario %in% c("Anull", "B3null")) %>%
  mutate(eval_type = forcats::fct_recode(eval_type,
                                         "Extraneous Factors" = "n_extra_0.7",
                                         "Mismatch Score" = "frob_n"),
         eval_type = factor(eval_type, levels = c("Extraneous Factors", "Mismatch Score")), # "RRMSE")),
         Method = forcats::fct_recode(method,
                                      "GFA-Rldsc" = "gfa-ldsc_pval_1_cn100:all",
                                      "GFA-Roracle" = "gfa-oracle_pval_1_cn100:all",
                                      "EBMF" = "gfa-none_pval_1_cn100:all",
                                      "MLFA" =  "ml-all_cn100:naf",
                                      "SPC" = "spc-5e-8:naf",
                                      "SVD" = "svd-0_5e-8:naf",
                                      "SVD-HT" = "svd-4_5e-8:naf",
                                      "Guide" = "guide-0_5e-8:naf"),
                                      #"FactorGo-1" = "factorgo-1_5e-8:naf"),
                                      #"FactorGo" = "factorgo-2_5e-8:naf"),
         Method = factor(Method, levels = c("GFA-Rldsc",  "GFA-Roracle",
                                            "EBMF", "MLFA",  #,  "FactorGO-1",
                                            "SPC", "SVD", "SVD-HT", "Guide"))) %>%
  ggplot(aes(x = scenario, y = eval,   color = Method), position = "dodge") +
  geom_boxplot(outlier.size = 1) +
  #geom_violin() +
  xlab("Scenario") +
  ylab("          Score           Num. Extra")+
  scale_color_manual( values = as.vector(palette.colors(palette = "Tableau 10", n = 10))[-5]) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 15)) +
  facet_wrap(~eval_type, ncol = 1, scale = "free_y",
  )

ggsave(p_null, file = paste0(dir, "/sims_null.png"), width = 10, height = 5, units = "in", dpi = 300)

## Blood Cell mimicking Scenario
# Effect of Condition Number on GFA
# Figure S10
p_bc_cn <- res_df_long %>%
  filter(method %in% c( "gfa-ldsc_pval_1_cn100:all", "gfa-oracle_pval_1_cn100:all",
                       "gfa-ldsc_pval_1_cn1000:all", "gfa-oracle_pval_1_cn1000:all",
                       "gfa-ldsc_pval_1_cn10000:all", "gfa-oracle_pval_1_cn10000:all"),
         eval_type %in% c("n_disc_0.9", "n_extra_0.7", "frob_n"),
         scenario %in% c("C")) %>%
  mutate(eval_type = forcats::fct_recode(eval_type,  "Discovered Factors" = "n_disc_0.9",
                                         "Extraneous Factors" = "n_extra_0.7",
                                         "Mismatch Score" = "frob_n"),
         eval_type = factor(eval_type, levels = c("Discovered Factors", "Extraneous Factors", "Mismatch Score")), # "RRMSE")),
         Method = forcats::fct_recode(method,
                                      "GFA-Rldsc_CN100" = "gfa-ldsc_pval_1_cn100:all",
                                      "GFA-Roracle_CN100" = "gfa-oracle_pval_1_cn100:all",
                                      #"EBMF" = "gfa-none_pval_1_cn100:all",
                                      "GFA-Rldsc_CN1000" = "gfa-ldsc_pval_1_cn1000:all",
                                      "GFA-Roracle_CN1000" = "gfa-oracle_pval_1_cn1000:all",
                                      "GFA-Rldsc_CN10000" = "gfa-ldsc_pval_1_cn10000:all",
                                      "GFA-Roracle_CN10000" = "gfa-oracle_pval_1_cn10000:all")) %>%
         # Method = factor(Method, levels = c("GFA_Rldsc",  "GFA_Roracle",
         #                                    "EBMF", "MLFA", "FactorGO-2", "FactorGO-1",
         #                                    "SPC", "SVD", "SVD-HT", "Guide"))) %>%
  ggplot(aes(x = scenario, y = eval,   color = Method), position = "dodge") +
  geom_boxplot() +
  xlab("Scenario") +
  ylab("          Score           Num. Extra       Num. Discovered")+
  scale_color_manual( values = as.vector(palette.colors(palette = "Tableau 10", n = 7))) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 15)) +
  facet_wrap(~eval_type, ncol = 1, scale = "free_y",
  )

ggsave(p_bc_cn, file = paste0(dir, "/sims_gfa_cn.png"), width = 8, height = 5, units = "in", dpi = 300)

## Big GFA Comparison
# Figure S6
p_gfa_compare <- res_df_long %>%
  filter(method %in% c( "gfa-ldsc_pval_1_cn100:all", "gfa-pt0.05_pval_1_cn100:all",
                        "gfa-ldsc_random_1_cn100:all", "gfa-pt0.05_random_1_cn100:all",
                        "gfa-ldsc_pval_1e-3_cn100:all", "gfa-pt0.05_pval_1e-3_cn100:all",
                        "gfa-ldsc_random_1e-3_cn100:all", "gfa-pt0.05_random_1e-3_cn100:all"),
         eval_type %in% c("n_disc_0.9", "n_extra_0.7", "frob_n"),
         scenario %in% c("A", "B1", "B2", "B3")) %>% #, "eval.rrmse")) %>%
  mutate(eval_type = forcats::fct_recode(eval_type,  "Discovered Factors" = "n_disc_0.9",
                                         "Extraneous Factors" = "n_extra_0.7",
                                         "Mismatch Score" = "frob_n"),
         eval_type = factor(eval_type, levels = c("Discovered Factors", "Extraneous Factors", "Mismatch Score")), # "RRMSE")),
         Method = forcats::fct_recode(method,
                                      "GFA-Rldsc-pval-1" = "gfa-ldsc_pval_1_cn100:all",
                                      "GFA-Rpt-pval-1" = "gfa-pt0.05_pval_1_cn100:all",
                                      "GFA-Rldsc-pval-1e-3" = "gfa-ldsc_pval_1e-3_cn100:all",
                                      "GFA-Rpt-pval-1e-3" =  "gfa-pt0.05_pval_1e-3_cn100:all",
                                      "GFA-Rldsc-random-1" = "gfa-ldsc_random_1_cn100:all",
                                      "GFA-Rpt-random-1" = "gfa-pt0.05_random_1_cn100:all",
                                      "GFA-Rldsc-random-1e-3" = "gfa-ldsc_random_1e-3_cn100:all",
                                      "GFA-Rpt-random-1e-3" =  "gfa-pt0.05_random_1e-3_cn100:all")) %>%
  ggplot(aes(x = scenario, y = eval,   color = Method), position = "dodge") +
  geom_boxplot() +
  xlab("Scenario") +
  ylab("          Score           Num. Extra       Num. Discovered")+
  scale_color_manual( values = as.vector(palette.colors(palette = "Tableau 10", n = 8))) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 15)) +
  facet_wrap(~eval_type, ncol = 1, scale = "free_y",
  )
ggsave(p_gfa_compare, file = paste0(dir, "/sims_compgfa.png"), width = 12, height = 5, units = "in", dpi = 300)


## Figures S11-S16
for(m in c("svd-0_5e-8:", "svd-4_5e-8:", "spc-5e-8:", "ml-all_cn100:", "guide-0_5e-8:", "factorgo-2_5e-8:")){
  lab <- case_match(m, "svd-0_5e-8:" ~ "SVD",
                    "svd-4_5e-8:" ~ "SVD-HT",
                    "spc-5e-8:" ~ "SPC",
                    "ml-all_cn100:" ~ "MLFA",
                    "guide-0_5e-8:" ~ "Guide",
                    "factorgo-2_5e-8:" ~ "FactorGo")
  p_nfactor  <- res_df_long %>%
                filter(stringr::str_starts(method, m)) %>%
                filter(eval_type %in% c("n_disc_0.9","n_extra_0.7", "frob_n"),
                      scenario %in% c("A", "B1", "B2", "B3", "C"),
                      ndense == 3 | scenario %in% c("Anull", "B3null", "C"),
                      method != "factorgo-2_5e-8:guess") %>%
                mutate(eval_type = forcats::fct_recode(eval_type,  "Discovered Factors" = "n_disc_0.9",
                                         "Extraneous Factors" = "n_extra_0.7",
                                         "Mismatch Score" = "frob_n"),
                       eval_type = factor(eval_type, levels = c("Discovered Factors", "Extraneous Factors", "Mismatch Score")),
                       num_factors = stringr::str_replace(method, m, ""),
                       num_factors = case_match(num_factors, "naf" ~ "Acceleration Factor",
                                  "nkaiser" ~ "Kaiser",
                                  "noc" ~ "Opt Coordinates",
                                  "six" ~ "Six",
                                  "twelve" ~ "Twelve")) %>%
                ggplot(aes(x = scenario, y = eval,   color = num_factors), position = "dodge") +
                geom_boxplot() +
                xlab("Scenario") +
                ylab("          Score           Num. Extra       Num. Discovered")+
                scale_color_manual( values = as.vector(palette.colors(palette = "Tableau 10", n = 9))) +
                theme_bw() +
                ggtitle(paste0("Choice of Factor Count for ", lab)) +
                theme(axis.text.y = element_text(size = 10),
                        axis.text.x = element_text(size = 20),
                        axis.title.x = element_text(size = 15),
                        axis.title.y = element_text(size = 12),
                        legend.text = element_text(size = 12),
                        legend.title = element_blank(),
                        strip.text = element_text(size = 15)) +
                facet_wrap(~eval_type, ncol = 1, scale = "free_y")
   ggsave(p_nfactor, file = paste0(dir, "/sims_nfactor_", lab, ".png"), width = 12, height = 5, units = "in", dpi = 300)
}


### Supplementary Plots Summarizing Simulations


## S1 Distribution of number of GW-sig variants for Scenario A3
files <- paste0("../../results/simulation_results/simdata/simdata.A.3.", 0:99, ".RDS")

nsig <- purrr::map(files, function(f){
  x <- readRDS(f)$nsig 
}) %>% unlist()
p <- ggplot(data.frame(n = nsig)) + geom_histogram(aes(x = n)) + xlab("Number of Independent Genome-Wide Significant Variants") + theme_bw() + theme(axis.text = element_text(size = 12))
ggsave(p, file = paste0(dir, "/nsig_dist_scenarioA3.png"), width = 6, height = 5, units = "in", dpi = 300)

## S2 Distribution of number of GW-sig variants for CAD/T2D traits
files <- paste0("../../data/gfa_intermediate_data/metab_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0_0.", 1:22, ".RDS")
t <- qnorm(1- ((5e-8)/2))
nsig <- purrr::map_dfr(files, function(f){
  x <- readRDS(f)
  Z <- select(x, ends_with(".z"))
  nsig <- apply(abs(Z), 2, function(xx){sum(xx > t)})
  return(nsig)
}) 
nsig <- colSums(nsig)
p <- ggplot(data.frame(n = nsig)) + geom_histogram(aes(x = n)) + xlab("Number of Independent Genome-Wide Significant Variants") + theme_bw()  + theme(axis.text = element_text(size = 12))
ggsave(p, file = paste0(dir, "/nsig_dist_metab.png"), width = 6, height = 5, units = "in", dpi = 300)


## Figure S4 Causal variant effect distribution
files <- paste0("../../results/simulation_results/simdata/simdata.A.3.", 0:99, ".RDS")

quant <- purrr::map2_dfr(files, 1:length(files), function(f, i){
  q <- readRDS(f)$quant_df %>% tidyr::pivot_longer(cols = paste0("x", 1:50))
  return(q)
}) %>% group_by(quantile) %>% summarize(med = median(value), q25 = quantile(value, prob  = 0.25), q75 = quantile(value, prob = 0.75))
p <- ggplot(quant) + geom_point(aes(x = quantile, y = med)) + 
  geom_errorbar(aes(ymin = q25, ymax = q75, x = quantile)) + 
  scale_y_log10() +
xlab("h2 percentile") + ylab("Effect Size (%h2)") +
theme_bw() +
theme(axis.text = element_text(size = 12), 
      axis.title = element_text(size = 15))
ggsave(p, file = paste0(dir, "/eff_dist_scenarioA3.png"), width = 6, height = 5, units = "in", dpi = 300)

## Figure S4 Genetic Correlation vs Residual Correlation

R <- readRDS("bc_R_estimate.R_ldsc.RDS")
Rg <- cov2cor(R$Rg)
Re <- cov2cor(R$R)
png(paste0(dir, "/bc_cg_v_ce.png"), width = 5.5, height = 500, res = 300, units = "in")
plot(Rg[lower.tri(Rg)], Re[lower.tri(Re)], xlab = "Pairwise Genetic Correlation", ylab = "Pairwise Residual Correlation")
abline(0, 1)
dev.off()


#### Figure S5 GW Sig variants in Scenario C vs in real BC data.
files <- paste0("../../data/gfa_intermediate_data/bc_zmat.ldpruned_r20.01_kb1000_pvalue_jitter0_0.", 1:22, ".RDS")
t <- qnorm(1- ((5e-8)/2))
nsig_bc_traits <- purrr::map_dfr(files, function(f){
  x <- readRDS(f)
  Z <- select(x, ends_with(".z"))
  nsig <- apply(abs(Z), 2, function(xx){sum(xx > t)})
  return(nsig)
}) 
df <- data.frame(trait = names(nsig_bc_traits), 
                 nsig_data = colSums(nsig_bc_traits))

files <- paste0("../../results/simulation_results/simdata/simdata.bc.3.", 0:99, ".RDS")
traits <- paste0(readRDS("bc_R_estimate.R_ldsc.RDS")$names, ".z")
nsig_sims <- purrr::map_dfr(files, function(f){
  x <- readRDS(f)$nsig 
  data.frame(trait = traits, value = x)
}) %>% group_by(trait) %>% summarize(med = median(value), q25 = quantile(value, prob  = 0.25), q75 = quantile(value, prob = 0.75))
df <- full_join(df, nsig_sims, by = "trait")
p <- ggplot(df) + geom_point(aes(x = nsig_data, y = med)) + 
  geom_errorbar(aes(ymin = q25, ymax = q75, x = nsig_data)) +
  geom_abline(slope = 1, intercept = 0) + 
  xlab("Number of Genome-Wide Significant Variants\nBlood Cell Data") + 
  ylab("Number of Genome-Wide Significant Variants\nSimulated Data") +
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15))
ggsave(p, file = paste0(dir, "/bc_ngwsig.png"), width = 6, height = 5, units = "in", dpi = 300)

