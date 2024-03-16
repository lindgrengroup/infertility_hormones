# Author: Samvida S. Venkatesh
# Date: 02/10/23

library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/selection_analyses"

INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis2",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5",
                   "male_infertility")

# Plot BetaScan2 scores vs -log10P windows of infertility and hormone loci ----

BETA_HIGH <- 4.1207

# Function to extract scores for all variants within 10 kb of lead variant 

WINDOW_SIZE <- 10000 # bp
getWindowScores <- function (beta_dat, chr_lead, pos_lead) {
  beta_window <- beta_dat %>% filter(CHR == paste0("chr", chr_lead) &
                                      POS >= pos_lead - WINDOW_SIZE &
                                      POS <= pos_lead + WINDOW_SIZE) %>%
    mutate(status = ifelse(BetaScan2_Std >= BETA_HIGH, "signif", "non_sig")) %>%
    rename(BP_pos = POS, score = BetaScan2_Std)
  return (beta_window)
}

# Function for p-value vs SDS plots ----

col_palette <- c("#A4A4A4", "#D35C79", "#009593", "#000000")
names(col_palette) <- c("non_sig", "signif_F", "signif_M", "lead")

scatterPlot <- function (locus_dat) {
  plot_dat <- locus_dat %>%
    mutate(status = ifelse(status == "signif" & sex_strata == "F", "signif_F",
                           ifelse(status == "signif" & sex_strata == "M", "signif_M",
                                  status)),
           neglogp = -log10(PVALUE))
  
  scatter_plot <- ggplot(plot_dat, aes(x = neglogp, y = score),
                           fill = status, colour = status) +
    geom_point(aes(fill = status, colour = status)) +
    geom_hline(yintercept = BETA_HIGH, linetype = "dashed") +
    scale_colour_manual(values = col_palette, guide = "none") +
    scale_fill_manual(values = col_palette, guide = "none") +
    labs(x = "-log10(PVALUE)", y = "BetaScan2") +
    theme(axis.text = element_text(size = 6, family = "Arial"),
          axis.title = element_blank()) 
  return (scatter_plot)
}

# Get plots for the four infertility loci of interest:

rsid <- "rs11692588"
chr_id <- 2
pos_id <- 11544358
tt <- "female_infertility_analysis5"
dat <- read.table(paste0(mainpath, "/data/", tt, "_EUR_betascan2_scores.txt"),
                  sep = " ", header = T, stringsAsFactors = F)
window_scores <- getWindowScores(dat, chr_id, pos_id)
window_scores$sex_strata <- "F"
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_beta_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)

rsid <- "rs10165819"
chr_id <- 2
pos_id <- 11581956
tt <- "female_infertility_analysis1"
dat <- read.table(paste0(mainpath, "/data/", tt, "_EUR_betascan2_scores.txt"),
                  sep = " ", header = T, stringsAsFactors = F)
window_scores <- getWindowScores(dat, chr_id, pos_id)
window_scores$sex_strata <- "F"
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_beta_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)

rsid <- "rs72827480"
chr_id <- 2
pos_id <- 120388925
tt <- "female_infertility_analysis3"
dat <- read.table(paste0(mainpath, "/data/", tt, "_EUR_betascan2_scores.txt"),
                  sep = " ", header = T, stringsAsFactors = F)
window_scores <- getWindowScores(dat, chr_id, pos_id)
window_scores$sex_strata <- "F"
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_beta_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)

rsid <- "rs150639836"
chr_id <- 10
pos_id <- 53879806
tt <- "male_infertility"
dat <- read.table(paste0(mainpath, "/data/", tt, "_EUR_betascan2_scores.txt"),
                  sep = " ", header = T, stringsAsFactors = F)
window_scores <- getWindowScores(dat, chr_id, pos_id)
window_scores$sex_strata <- "M"
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_beta_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)
