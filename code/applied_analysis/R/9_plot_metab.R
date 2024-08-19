library(dplyr)
library(ggplot2)
library(stringr)
library(GFA)
library(viridis)


sessionInfo()




plotf <- function(fit, o, fct_order, flip_sign = FALSE){
  if(missing(o)){ o <- seq_along(fit$names)}
  if(missing(fct_order)) fct_order <- seq_along(ncol(fit1$F_hat))
  if(flip_sign){
    myF <- lapply(seq(ncol(fit$F_hat)), function(i){
      ii <- which.max(fit$gfa_pve$pve[,i])
      sign(fit$F_hat[ii,i])*fit$F_hat[,i]}) %>% do.call(cbind, .)
  }else{
    myF <- fit$F_hat
  }
  plot_factors(myF[o,fct_order], row_names = fit$names[o]) +
    scale_fill_gradient2(low = viridis(4)[3], mid = "white", high = viridis(4)[1],
                         name = "Factor Effect") +
    xlab("Factor") + ylab("Trait") +
    theme(axis.text.x  = element_text(size=14),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size=15),
          strip.text = element_text(size=15))
}

gfa_fit <- readRDS(snakemake@input[["gfa_fit"]])

####
# convert to liability scale for binary traits
N <- gfa_fit$scale^2 # originally entered N
binary_traits <- c( "van-der-Harst_2017_29212778__CAD",
                    "Xue_2018_30054458__type-2-diabetes",
                    "Liu_2019_30643251__smoking-initiation",
                    "Elsworth_2018_0__Sleep-apnoea",
                    "Elsworth_2018_0__Loneliness-isolation",
                    "Trubetskoy_2022_35396580__Schizophrenia",
                    "Demontis_2019_30478444__ADHD",
                    "Mullins_2021_34002096__Bipolar-disorder",
                    "Howard_2019_30718901__major-depression")

N_case <- pop_prev <- rep(NA, length(N))
N_case[match(binary_traits, gfa_fit$names)] <- c(122733, 62892, 311629,
                                               2320, 82436, 76755, 20183,
                                               41917, 170756)
pop_prev[match(binary_traits, gfa_fit$names)] <- c(0.11, 0.1, 0.13, 0.1, 0.18, 0.005,
                                                 0.08, 0.044, 0.083)

convert_to_liab_scale <- function(gfa_res, N, N_case, pop_prev){
  new_scale <- gfa_res$scale
  ix <- which(!is.na(N_case))
  if(length(ix) > 0) new_scale[ix] <- GFA:::binary_const(N = N[ix], N_case = N_case[ix], pop_prev = pop_prev[ix])
  gfa_res_new <- gfa_wrapup(fit = gfa_res$fit, method = gfa_res$method, scale = new_scale)
  gfa_res_new$names <- gfa_res$names
  return(gfa_res_new)
}

gfa_fit_liab <- convert_to_liab_scale(gfa_fit, N, N_case, pop_prev)
##########


short_name <- str_split(gfa_fit$names, "__") %>% purrr::map(2) %>% unlist()
name_dict <- data.frame(orig_name = gfa_fit$name,
                        name = short_name) %>%
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
gfa_fit$names <- gfa_fit_liab$names <- name_dict$print_name[match(gfa_fit$names, name_dict$orig_name)]
print_names <- c( "CAD","T2D",    "BMI" ,   "BFP" ,
                 "Waist Circ.", "Hip Circ.", "Body Size Age 10",
                 "BMR", "Height",           "LDL" ,             "HDL",              "Trig",
                  "SBP"       ,       "DBP"     ,         "CRP" ,             "Smoking",
                  "Alc/Wk"   ,        "Insomnia"  ,       "Sleep Apnoea",     "Loneliness"    ,
                   "ADHD"     ,        "MDD"   ,           "BPD"   ,           "Scz"      )
o <- match(print_names, name_dict$print_name)
fct_order <- c(1, 5, 16, 12, 2,  4, 7, 10, 3, 14, 6, 9, 11, 8, 13, 15)

p1_liab <- plotf(gfa_fit_liab, o, fct_order = fct_order, flip_sign = TRUE)

percents <- round(100*gfa_fit$gfa_pve$pve[o, fct_order], digits = 1) %>% as.vector()
value <- as.vector(gfa_fit$F_hat[o, fct_order])
pdf <- data.frame(label = paste0(percents, "%"),
                  percent = percents,
                  value = value,
                  x = rep(1:16, each =24),
                  y = rep(c(1:24), 16)) %>%
        mutate(color = case_when(abs(value) < 0.6 ~ "black",
                                 TRUE ~ "white")) %>% filter(percents >= 0.5)

p2_liab <- p1_liab + geom_text(data = pdf ,
                                 aes(x = x, y = y, label = label), color = pdf$color, size = 2.8)
ggsave(p2_liab, file = snakemake@output[["out"]],create.dir = TRUE, height = 6, width = 12, units = "in", dpi = 300)

