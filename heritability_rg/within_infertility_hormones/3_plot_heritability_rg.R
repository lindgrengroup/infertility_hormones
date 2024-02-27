# Author: Samvida S. Venkatesh
# Date: 16/03/2023

library(tidyverse)
theme_set(theme_bw())

# Read data ----

HORMONES <- c("FSH_F", "FSH_M",
              "LH_F", "LH_M",
              "Oestradiol_F", "Oestradiol_M", 
              "Testosterone_F", "Testosterone_M")
# Remove male infertility for all-anc
INFERTILITY <- c("female_infertility_analysis1",
                 "female_infertility_analysis2",
                 "female_infertility_analysis3",
                 "female_infertility_analysis4",
                 "female_infertility_analysis5",
                 "male_infertility")

EUR_rg <- read.table("EUR_genetic_correlations.txt",
                     sep = "\t", header = T, stringsAsFactors = F)

EUR_h2 <- read.table("EUR_heritability.txt",
                     sep = "\t", header = T, stringsAsFactors = F)

# Add a set of columns for correlations with themselves
self_corr <- data.frame(trait_1 = c(HORMONES, INFERTILITY),
                        trait_2 = c(HORMONES, INFERTILITY),
                        rg = 1,
                        se = 0,
                        pval = 0)
EUR_rg <- bind_rows(EUR_rg, self_corr)

# Plot heritability ----

custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

getHeritPlot <- function (herit_df) {
  h2_plot <- ggplot(herit_df, 
                    aes(x = heritability, y = strata,
                        colour = sex_strata, fill = sex_strata)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        colour = sex_strata, fill = sex_strata),
                    size = 0.7,
                    position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = custom_three_diverge, guide = "none") +
    scale_colour_manual(values = custom_three_diverge, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme(axis.title.y = element_blank())
  return (h2_plot)
}

EUR_h2 <- EUR_h2 %>%
  mutate(lci = heritability - 1.96*se,
         uci = heritability + 1.96*se,
         strata = factor(strata, levels = c(HORMONES, INFERTILITY)))

hormones_EUR_h2 <- EUR_h2 %>%
  filter(strata %in% HORMONES) %>%
  mutate(sex_strata = gsub(".*_", "", strata))
hormones_EUR_h2$sex_strata[hormones_EUR_h2$sex_strata == "comb"] <- "sex_comb"

ggsave("plots/hormones_EUR_h2.png", getHeritPlot(hormones_EUR_h2), 
       height = 7, width = 7, units = "in")

infertility_EUR_h2 <- EUR_h2 %>%
  filter(strata %in% INFERTILITY) %>%
  mutate(sex_strata = ifelse(grepl("^female", strata), "F", "M"))

ggsave("plots/infertility_EUR_h2.png", getHeritPlot(infertility_EUR_h2), 
       height = 7, width = 7, units = "in")

# Plot correlations ----

alpha_vals <- c(0.5, 1)
names(alpha_vals) <- c("no", "signif")

getCorrPlot <- function (corr_df_long) {
  plot_rg <- ggplot(corr_df_long,
                             aes(x = trait_1, y = trait_2, fill = rg)) +
    geom_tile(aes(alpha = signif), colour = "black") +
    geom_text(data = corr_df_long %>% filter(signif == "signif"),
              aes(label = rg_lab), colour = "black") +
    scale_fill_gradient2(low = "#009593", mid = "#FFFFCC", high = "#D35C79", 
                         limits = c(-1, 1), midpoint = 0) +
    scale_alpha_manual(values = alpha_vals, guide = "none") +
    coord_fixed() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 30,
                                     hjust = 1))
  return (plot_rg)
}

# Hormones ---

hormones_EUR_rg <- EUR_rg %>% 
  filter(trait_1 %in% HORMONES &
           trait_2 %in% HORMONES & 
           !is.na(rg) & rg <= 1 & rg >= -1) %>%
  mutate(signif = ifelse(pval < 0.05, "signif", "no"),
         trait_1 = factor(trait_1, levels = HORMONES),
         trait_2 = factor(trait_2, levels = rev(HORMONES)),
         rg_lab = signif(rg, 2))

ggsave("plots/hormones_EUR_rg.png", getCorrPlot(hormones_EUR_rg), 
       height = 7, width = 7, units = "in")

# Infertility ----

inf_EUR_rg <- EUR_rg %>% 
  filter(trait_1 %in% INFERTILITY &
           trait_2 %in% INFERTILITY &
           trait_1 != "male_infertility" & trait_2 != "male_infertility" &
           !is.na(rg) & rg <= 1 & rg >= -1) %>%
  mutate(signif = ifelse(pval < 0.05, "signif", "no"),
         trait_1 = factor(trait_1, levels = INFERTILITY),
         trait_2 = factor(trait_2, levels = rev(INFERTILITY)),
         rg_lab = signif(rg, 2))

ggsave("plots/infertility_EUR_rg.png", getCorrPlot(inf_EUR_rg), 
       height = 7, width = 7, units = "in")

# Hormones vs infertility ----

hinf_EUR_rg <- EUR_rg %>% 
  filter(trait_1 %in% HORMONES & 
           trait_2 %in% INFERTILITY & 
           !is.na(rg) & rg <= 1 & rg >= -1) %>%
  mutate(signif = ifelse(pval < 0.05, "signif", "no"),
         trait_1 = factor(trait_1, levels = HORMONES),
         trait_2 = factor(trait_2, levels = INFERTILITY),
         rg_lab = signif(rg, 2))

ggsave("plots/hormones_vs_inf_EUR_rg.png", getCorrPlot(hinf_EUR_rg), 
       height = 7, width = 7, units = "in")

getCorrPlotPoster <- function (corr_df_long) {
  plot_rg <- ggplot(corr_df_long,
                    aes(x = trait_1, y = trait_2, fill = rg)) +
    geom_tile(aes(alpha = signif), colour = "black") +
    scale_fill_gradient2(low = "#009593", mid = "#FFFFCC", high = "#D35C79", 
                         na.value = "grey",
                         limits = c(-1, 1), midpoint = 0) +
    scale_alpha_manual(values = alpha_vals, guide = "none") +
    coord_fixed() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  return (plot_rg)
}


for_poster <- EUR_rg %>% 
  filter(trait_1 %in% c("FSH_F", "LH_F", "Oestradiol_F", "Progesterone_F", "Testosterone_F") & 
           trait_2 %in% c("female_infertility_analysis1", 
                          "female_infertility_analysis3",
                          "female_infertility_analysis5")) %>%
  mutate(signif = ifelse(pval < 0.05, "signif", "no"),
         trait_1 = factor(trait_1, levels = c("FSH_F", "LH_F", "Oestradiol_F",
                                              "Progesterone_F", "Testosterone_F")),
         trait_2 = factor(trait_2, levels = rev(INFERTILITY)),
         pval_lab = paste0("P = ", signif(pval, 2)))

ggsave("plots/hormones_vs_inf_for_poster.tiff", getCorrPlotPoster(for_poster), 
       height = 9, width = 22, units = "cm")


