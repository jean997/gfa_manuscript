library(dplyr)
library(ggplot2)
library(stringr)
library(GFA)
library(viridis)


sessionInfo()

source("R/order_scale_factors.R")

res <- readRDS(snakemake@input[["gfa_fit"]])

##########


names_simple <- str_split(res$names, "__") %>% purrr::map(2) %>% unlist()
name_dict <- data.frame(orig_name = res$name,
                        name = names_simple) %>%
              mutate(print_name = case_when(name == "type-2-diabetes" ~ "T2D",
                                            name == "smoking-initiation" ~ "Smoking",
                                            name == "Alcoholic-drinks-per-week" ~ "Alc/Wk",
                                            name == "Sleeplessness-insomnia" ~ "Insomnia",
                                            name == "Sleep-apnoea" ~ "Sleep Apnoea",
                                            name == "Waist-circumference"  ~ "Waist Circ.",
                                            name == "Body-fat-percentage" ~ "BFP",
                                            name == "LDL-cholesterol" ~ "LDL",
                                            name == "HDL-cholesterol" ~ "HDL",
                                            name == "triglycerides" ~ "Trig",
                                            name == "Hip-circumference" ~ "Hip Circ.",
                                            name == "Comparative-body-size-at-age-10" ~ "Body Size Age 10",
                                            name == "Loneliness-isolation" ~ "Loneliness",
                                            name == "Schizophrenia"   ~ "Scz",
                                            name == "Bipolar-disorder" ~ "BPD",
                                            name == "major-depression"  ~ "MDD",
                                            TRUE ~ name
                                            ))
res$names <-  name_dict$print_name[match(res$names, name_dict$orig_name)]
trait_order <- c( "CAD","T2D",    "BMI" ,   "BFP" ,
                 "Waist Circ.", "Hip Circ.", "Body Size Age 10",
                 "BMR", "Height",           "LDL" ,             "HDL",              "Trig",
                  "SBP"       ,       "DBP"     ,         "CRP" ,             "Smoking",
                  "Alc/Wk"   ,        "Insomnia"  ,       "Sleep Apnoea",     "Loneliness"    ,
                   "ADHD"     ,        "MDD"   ,           "BPD"   ,           "Scz"      )
o <- match(trait_order, res$names)
nf <- ncol(res$F_hat)

F_hat <- order_and_scale(res$F_hat[o,])
fct_order <- F_hat$order

# put cad and t2d at the end
fct_order <- fct_order[c(3:nf, 1, 2)]
F_hat$F_hat <- F_hat$F_hat[,c(3:nf, 1, 2)]


p1 <- F_hat$F_hat %>%
  plot_factors(row_names = res$names[o]) +
  scale_fill_gradient2(low = viridis(4)[3], mid = "white", high = viridis(4)[1],
                       name = "Effect") +
  xlab("Factor") + ylab("Trait") +
  theme(axis.text.x  = element_text(size=14, angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size=15),
        strip.text = element_text(size=15))

percents <- round(100*res$gfa_pve$pve[o, fct_order], digits = 1) %>% as.vector()
value <- as.vector(F_hat$F_hat)
pdf <- data.frame(label = paste0(percents, "%"),
                  percent = percents,
                  value = value,
                  x = rep(1:nf, each =24),
                  y = rep(c(1:24), nf)) %>%
        mutate(color = case_when(abs(value) < 0.6 ~ "black",
                                 TRUE ~ "white")) %>% filter(percents >= 0.5)

p2 <- p1 + geom_text(data = pdf ,
                     aes(x = x, y = y, label = label), color = pdf$color, size = 2.8)
ggsave(p2, file = snakemake@output[["out"]],create.dir = TRUE, height = 6, width = 12, units = "in", dpi = 300)

## compute pve credible intervals
nsamp <- as.numeric(snakemake@params[["nsamp"]])
level <- as.numeric(snakemake@params[["level"]])
pve_int <- GFA:::gfa_credints(res, nsamp = nsamp, level = level, type = "pve")

# reorder to match publication
pve_int$pve_lower <- pve_int$pve_lower[o, fct_order]
pve_int$pve_upper <- pve_int$pve_upper[o, fct_order]
pve_int$pve_est <-  res$gfa_pve$pve[o, fct_order]
pve_int$names <- res$names[o]
pve_int$F_hat <- F_hat$F_hat

saveRDS(pve_int, file = snakemake@output[["pve_ci"]])
