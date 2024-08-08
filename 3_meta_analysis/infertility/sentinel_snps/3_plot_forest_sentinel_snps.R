# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(tidyverse)
library(ggpubr)
library(scales)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility"

F_INFERT <- c("female_infertility_analysis1", "female_infertility_analysis2",
              "female_infertility_analysis3", "female_infertility_analysis4",
              "female_infertility_analysis5")

ANC_GPS <- c("eur", "all")

# Read data on sentinel SNPs so we know which SNP came from which analysis ----

sentinel_snps <- lapply(F_INFERT, function (finf) {
  if (finf == "female_infertility_analysis1") {
    f_snps <- read.table(paste0(mainpath, 
                                "/infertility_meta_mvp/lead_snps/female_infertility_analysis1_all_sentinel_SNPs.txt"),
                         sep = "\t", header = T, stringsAsFactors = F)$ID
  } else {
    f_snps <- read.table(paste0("/well/lindgren/laura/projects/infertility/sentinel_snps/", 
                                finf, "_all_sentinel_SNPs.txt"), 
                         sep = "\t", header = T, stringsAsFactors = F)$MarkerName
  }
  return (f_snps)
})
names(sentinel_snps) <- F_INFERT

# Wrangle data ----

# Standardise columns to plot
to_numeric <- c("Freq1", "FreqSE", "Effect", "StdErr", "P.value",
                "HetISq", "HetChiSq", "HetDf", "HetPVal",
                "GENPOS", "MAF", "BETA", "SE", "PVALUE")

f_dat <- lapply(F_INFERT, function (finf) {
  per_anc <- lapply(ANC_GPS, function (anc) {
    # print(paste0(finf, "_", anc))
    df <- read.table(paste0(mainpath, "/infertility_sumstats/all_lead_snps_", 
                            finf, "_", anc, ".txt"),
                     sep = "\t", header = T, stringsAsFactors = F)
    df <- df %>%
      mutate(across(any_of(to_numeric), as.numeric),
             pheno = finf, ancestry = anc)
    
    if ("MarkerName" %in% colnames(df)) {
      df <- df %>%
        rename(ID = MarkerName, BETA = Effect, SE = StdErr, PVALUE = P.value) 
    }
    df <- df %>%
      mutate(flipped_BETA = ifelse(Freq1 > 0.5, -BETA, BETA)) %>%
      select(all_of(c("ID", "pheno", "ancestry", "flipped_BETA", "SE", "PVALUE")))
    return (df)
  })
  per_anc <- bind_rows(per_anc)
  return (per_anc)
})
f_dat <- bind_rows(f_dat)

f_dat <- f_dat %>%
  mutate(OR = exp(flipped_BETA), 
         LCI = exp(flipped_BETA-1.96*SE), 
         UCI = exp(flipped_BETA+1.96*SE),
         pheno = factor(pheno, levels = rev(F_INFERT)),
         gws = ifelse(PVALUE <= 5E-08, "yes", "no"),
         nominal_sig = ifelse(PVALUE <= 0.05/22, "yes", "no"))

# Forest plots ----

line_palette <- c("solid", "dashed")
names(line_palette) <- c("yes", "no")

f_plot_per_id <- lapply(unique(f_dat$ID), function (id_plot) {
  dat_plot <- f_dat %>% filter(ID == id_plot)
  fplot <- ggplot(dat_plot, aes(x = OR, y = pheno,
                                group = ancestry)) +
    geom_point(aes(shape = ancestry), 
               position = position_dodge(width = 0.7), 
               colour = "#D35C79", size = 1) +
    geom_errorbar(aes(xmin = LCI, xmax = UCI, linetype = gws), 
                  position = position_dodge(width = 0.7), 
                  colour = "#D35C79", size = 0.2) +
    scale_shape(guide = "none") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(trans = "log",
                       breaks = pretty_breaks(n = 3)) +
    scale_linetype_manual(values = line_palette, guide = "none") +
    scale_y_discrete(drop = F) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_blank())
  return (fplot)
})
names(f_plot_per_id) <- unique(f_dat$ID)

p1 <- ggarrange(plotlist = f_plot_per_id[sentinel_snps[[1]]], 
                nrow = 2, ncol = 4)
ggsave(paste0(mainpath, "/infertility_meta_plots/forest_plot_female_infertility_set1.png"), p1,
       width = 14, height = 6, units = "cm")
p2 <- ggarrange(plotlist = f_plot_per_id[sentinel_snps[[2]]], 
                nrow = 2, ncol = 4)
ggsave(paste0(mainpath, "/infertility_meta_plots/forest_plot_female_infertility_set2.png"), p2,
       width = 14, height = 6, units = "cm")
p3 <- ggarrange(plotlist = f_plot_per_id[sentinel_snps[[3]]], 
                nrow = 2, ncol = 4)
ggsave(paste0(mainpath, "/infertility_meta_plots/forest_plot_female_infertility_set3.png"), p3,
       width = 14, height = 6, units = "cm")
p4 <- ggarrange(plotlist = f_plot_per_id[sentinel_snps[[4]]], 
                nrow = 2, ncol = 4)
ggsave(paste0(mainpath, "/infertility_meta_plots/forest_plot_female_infertility_set4.png"), p4,
       width = 14, height = 6, units = "cm")
p5 <- ggarrange(plotlist = f_plot_per_id[sentinel_snps[[5]]], 
                nrow = 2, ncol = 4)
ggsave(paste0(mainpath, "/infertility_meta_plots/forest_plot_female_infertility_set5.png"), p5,
       width = 14, height = 6, units = "cm")
