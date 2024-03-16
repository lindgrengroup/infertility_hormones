# Author: Samvida S. Venkatesh
# Date: 02/10/23

library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/selection_analyses"
HORMONE_STRATA <- c("FSH_F", "FSH_M",
                    "LH_F", "LH_M",
                    "Oestradiol_F", "Oestradiol_M",
                    "Testosterone_F", "Testosterone_M")

INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis2",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5",
                   "male_infertility")

# Function to read and flip selection scores ----

getTSDS <- function (trait) {
  dat <- read.table(paste0(mainpath, "/data/", trait, "_EUR_sds_scores.txt"),
                    sep = " ", header = T, stringsAsFactors = F)
  
  # Align to trait positive alleles (i.e. hormone-increasing or infertility-risk increasing)
  tmp <- dat %>%
    mutate(flip_beta = BETA < 0,
           new_beta = ifelse(flip_beta, -BETA, BETA),
           new_tested_allele_freq = ifelse(flip_beta, 1 - Tested_Allele_Freq, Tested_Allele_Freq),
           new_tested_allele = ifelse(flip_beta, Other_Allele, Tested_Allele),
           new_other_allele = ifelse(flip_beta, Tested_Allele, Other_Allele))
  # Align tSDS to the new tested allele
  tmp <- tmp %>%
    mutate(derived_test_match = Derived_Allele == new_tested_allele,
           tSDS = ifelse(derived_test_match, SDS, -SDS))
  
  res <- tmp[, c("CHR", "POS", "RSID", 
                 "Ancestral_Allele", "Derived_Allele", "Derived_Allele_Freq",
                 "SDS", "tSDS", "PVALUE")]
  res$Tested_Allele_Freq <- tmp$new_tested_allele_freq
  res$Tested_Allele <- tmp$new_tested_allele
  res$Other_Allele <- tmp$new_other_allele
  res$BETA <- tmp$new_beta
  return (res)
}

# Plot rank-based mean tSDS ----

NSNPS_GP <- 1000
getTSDSPlot <- function (tsds_dat, trait) {
  # Rank by pval 
  tsds_dat <- tsds_dat %>% arrange(desc(PVALUE)) %>%
    mutate(pval_rank = 1:nrow(tsds_dat),
           pval_rank_gp = pval_rank %/% NSNPS_GP)
  
  for_plot <- tsds_dat %>% 
    group_by(pval_rank_gp) %>%
    mutate(mean_tsds = mean(tSDS, na.rm = T))
  
  tSDS_plot <- ggplot(for_plot, aes(x = pval_rank_gp,
                                    y = mean_tsds)) +
    geom_point(colour = "grey") +
    geom_hline(yintercept = 0) +
    labs(x = paste0(trait, " GWAS significance rank"),
         y = paste0(trait, "-increasing tSDS (bin average)"))
  return (tSDS_plot)
}

# Get mean tSDS for all trait-increasing alleles ----

GWS_THRESH <- 5E-08
mean_tsds <- lapply(c(HORMONE_STRATA, INFERT_STRATA), function (tt) {
  tsds_scores <- getTSDS(tt)
  # Plot
  ggsave(paste0(mainpath, "/plots/tSDS_", tt, ".png"), 
         getTSDSPlot(tsds_scores, tt))
  # Get mean
  gws_tsds_scores <- tsds_scores %>% filter(PVALUE <= GWS_THRESH)
  res <- data.frame(trait = tt,
                    mean_tsds = mean(gws_tsds_scores$tSDS),
                    sd_tsds = sd(gws_tsds_scores$tSDS))
  return (res)
})
mean_tsds <- bind_rows(mean_tsds)
write.table(mean_tsds, paste0(mainpath, "/results/trait_sds_scores.txt"),
            sep = "\t", row.names = F, quote = F)

# Plot tSDS scores in windows of infertility and hormone loci ----

SDS_HIGH <- 1.99441
SDS_LOW <- -1.97067

# Function to extract scores for all variants within 10 kb of lead variant 

WINDOW_SIZE <- 10000 # bp
getWindowScores <- function (tsds_dat, chr_lead, pos_lead) {
  sds_window <- tsds_dat %>% filter(CHR == paste0("chr", chr_lead) &
                                      POS >= pos_lead - WINDOW_SIZE &
                                      POS <= pos_lead + WINDOW_SIZE) %>%
    mutate(status = ifelse(tSDS >= SDS_HIGH | tSDS <= SDS_LOW, "signif", "non_sig")) %>%
    rename(BP_pos = POS, score = tSDS)
  return (sds_window)
}

# Function for mini-Manhattan plots within each lead variant locus ----

col_palette <- c("#A4A4A4", "#D35C79", "#009593", "#000000")
names(col_palette) <- c("non_sig", "signif_F", "signif_M", "lead")

locusPlot <- function (locus_dat, locus_name, 
                       locus_position) {
  plot_dat <- locus_dat %>%
    mutate(status = ifelse(status == "signif" & sex_strata == "F", "signif_F",
                           ifelse(status == "signif" & sex_strata == "M", "signif_M",
                                  status)))
  # Add text for lead variant
  if (locus_position %in% plot_dat$BP_pos) {
    for_lab <- plot_dat[plot_dat$BP_pos == locus_position, ]
    for_lab$label_id <- locus_name
    for_lab$status <- "lead"
  } else {
    for_lab <- data.frame(BP_pos = locus_position,
                          score = 0, status = "lead",
                          label_id = locus_name)
  }
  plot_dat <- plot_dat %>% filter(BP_pos != locus_position)
  plot_dat <- bind_rows(plot_dat, for_lab)
  
  minpos <- min(plot_dat$BP_pos)
  maxpos <- max(plot_dat$BP_pos)
  
  manhattan_plot <- ggplot(plot_dat, aes(x = BP_pos, y = score),
                           fill = status, colour = status) +
    geom_point(aes(fill = status, colour = status)) +
    geom_text_repel(data = for_lab, aes(label = label_id),
                    size = 2.5, family = "Arial") +
    geom_hline(yintercept = SDS_LOW, linetype = "dashed") +
    geom_hline(yintercept = SDS_HIGH, linetype = "dashed") +
    scale_colour_manual(values = col_palette, guide = "none") +
    scale_fill_manual(values = col_palette, guide = "none") +
    scale_x_continuous(limits = c(minpos, maxpos), breaks = pretty_breaks()) +
    labs(x = "Position (hg38)", y = "trait-SDS") +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(size = 6, family = "Arial"),
          axis.title = element_blank()) 
  return (manhattan_plot)
}

# Function for p-value vs SDS plots ----

scatterPlot <- function (locus_dat) {
  plot_dat <- locus_dat %>%
    mutate(status = ifelse(status == "signif" & sex_strata == "F", "signif_F",
                           ifelse(status == "signif" & sex_strata == "M", "signif_M",
                                  status)),
           neglogp = -log10(PVALUE))
  
  scatter_plot <- ggplot(plot_dat, aes(x = neglogp, y = score),
                           fill = status, colour = status) +
    geom_point(aes(fill = status, colour = status)) +
    geom_hline(yintercept = SDS_LOW, linetype = "dashed") +
    geom_hline(yintercept = SDS_HIGH, linetype = "dashed") +
    scale_colour_manual(values = col_palette, guide = "none") +
    scale_fill_manual(values = col_palette, guide = "none") +
    labs(x = "-log10(PVALUE)", y = "trait-SDS") +
    theme(axis.text = element_text(size = 6, family = "Arial"),
          axis.title = element_blank()) 
  return (scatter_plot)
}

# Get plots for the two infertility loci:
rsid <- "rs1964514"
chr_id <- 8
pos_id <- 109463457
tt <- "female_infertility_analysis1"
tsds_scores <- getTSDS(tt)
window_scores <- getWindowScores(tsds_scores, chr_id, pos_id)
window_scores$sex_strata <- "F"
lplot <- locusPlot(locus_dat = window_scores, 
                   locus_name = rsid, 
                   locus_position = pos_id)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds.png"),
       lplot,
       units = "cm", height = 4.5, width = 7)
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)


rsid <- "rs72676844"
chr_id <- 8
pos_id <- 108958927
tt <- "female_infertility_analysis4"
tsds_scores <- getTSDS(tt)
window_scores <- getWindowScores(tsds_scores, chr_id, pos_id)
window_scores$sex_strata <- "F"
lplot <- locusPlot(locus_dat = window_scores, 
                   locus_name = rsid, 
                   locus_position = pos_id)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds.png"),
       lplot,
       units = "cm", height = 4.5, width = 7)
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)

# Get plots for the four testosterone loci:
rsid <- "rs7578292"
chr_id <- 2
pos_id <- 136199163
tt <- "Testosterone_F"
tsds_scores <- getTSDS(tt)
window_scores <- getWindowScores(tsds_scores, chr_id, pos_id)
window_scores$sex_strata <- "F"
lplot <- locusPlot(locus_dat = window_scores, 
                   locus_name = rsid, 
                   locus_position = pos_id)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds.png"),
       lplot,
       units = "cm", height = 4.5, width = 7)
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)


rsid <- "rs1185977"
chr_id <- 6
pos_id <- 25828826 
tt <- "Testosterone_F"
tsds_scores <- getTSDS(tt)
window_scores <- getWindowScores(tsds_scores, chr_id, pos_id)
window_scores$sex_strata <- "F"
lplot <- locusPlot(locus_dat = window_scores, 
                   locus_name = rsid, 
                   locus_position = pos_id)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds.png"),
       lplot,
       units = "cm", height = 4.5, width = 7)
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)


rsid <- "rs112881196"
chr_id <- 2
pos_id <- 31757742  
tt <- "Testosterone_M"
tsds_scores <- getTSDS(tt)
window_scores <- getWindowScores(tsds_scores, chr_id, pos_id)
window_scores$sex_strata <- "M"
lplot <- locusPlot(locus_dat = window_scores, 
                   locus_name = rsid, 
                   locus_position = pos_id)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds.png"),
       lplot,
       units = "cm", height = 4.5, width = 7)
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)


rsid <- "rs8016626"
chr_id <- 14
pos_id <- 73641646  
tt <- "Testosterone_M"
tsds_scores <- getTSDS(tt)
window_scores <- getWindowScores(tsds_scores, chr_id, pos_id)
window_scores$sex_strata <- "M"
lplot <- locusPlot(locus_dat = window_scores, 
                   locus_name = rsid, 
                   locus_position = pos_id)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds.png"),
       lplot,
       units = "cm", height = 4.5, width = 7)
splot <- scatterPlot(locus_dat = window_scores)
ggsave(paste0(mainpath, "/plots/", rsid, "_tsds_scatter.png"),
       splot,
       units = "cm", height = 4.5, width = 4.5)


