# Author: Samvida S. Venkatesh
# Date: 16/03/2023

library(tidyverse)
library(RColorBrewer)
library(ggrepel)
theme_set(theme_bw())

# Read data ----

main_filepath <- "/well/lindgren/samvida/hormones_infertility/rg_hormones_infertility"

# Remove Oestradiol-M for all-anc
HORMONES <- c("FSH_F", "FSH_sex_comb",
              "LH_F", "LH_sex_comb",
              "Oestradiol_F", "Oestradiol_M", "Oestradiol_sex_comb",
              "Progesterone_F", "Progesterone_sex_comb",
              "Testosterone_F", "Testosterone_M", "Testosterone_sex_comb")
# Remove male infertility for all-anc
INFERTILITY <- c("female_infertility_analysis1",
                 "female_infertility_analysis2",
                 "female_infertility_analysis3",
                 "female_infertility_analysis4",
                 "female_infertility_analysis5",
                 "male_infertility")

# Categories for visualisation
celltypes_tissue <- read.table("/well/lindgren/samvida/hormones_infertility/ovary_celltype_markers_melody/celltype_categories.txt",
                               sep = "\t", header = T, stringsAsFactors = F)
celltypes_tissue <- celltypes_tissue %>%
  mutate(category = factor(category,
                                  levels = c("Endothelial", "Epithelial",
                                             "Granulosa", "Immune",
                                             "Smooth_Muscle", "Stroma_Theca",
                                             "Other")))

# European results
cts_EUR <- lapply(c(HORMONES, INFERTILITY), function (trait) {
  checkname <- paste0(main_filepath, "/ovary_celltypes_heritability/",
                      trait, "_EUR_ovary.cell_type_results.txt")
  if (file.exists(checkname)) {
    dat <- read.table(checkname,
                      sep = "\t", header = T)
    dat$trait <- trait
    dat$FDR_pval <- p.adjust(dat$Coefficient_P_value, method = "fdr")
  } else dat <- NULL
  return (dat)
})
cts_EUR <- bind_rows(cts_EUR)

# All ancestry results
# cts_all_anc <- lapply(c(HORMONES, INFERTILITY), function (trait) {
#   checkname <- paste0(main_filepath, "/ovary_celltypes_heritability/",
#                       trait, "_all_anc_ovary.cell_type_results.txt")
#   if (file.exists(checkname)) {
#     dat <- read.table(checkname,
#                       sep = "\t", header = T)
#     dat$trait <- trait
#     dat$FDR_pval <- p.adjust(dat$Coefficient_P_value, method = "fdr")
#   } else dat <- NULL
#   return (dat)
# })
# cts_all_anc <- bind_rows(cts_all_anc)

# Plot enrichment plots ----

PTHRESH <- 1.6 # This is -log10 for FDR < 5%
alpha_vals <- c(0.5, 1)
names(alpha_vals) <- c("no", "signif")

# Expand colour palette
colpal_use <- colorRampPalette(brewer.pal(9, "Set1"))(9)

getEnrichmentPlot <- function (df) {
  to_label <- df %>% filter(signif == "signif") %>%
    group_by(category, trait) %>%
    slice_max(order_by = neglogp, n = 2)
  
  res_plot <- ggplot(df, aes(x = Name, y = neglogp)) +
    #facet_wrap(~trait, scales = "free_y") +
    geom_point(aes(size = neglogp, alpha = signif,
                   fill = category, colour = category)) +
    geom_text_repel(data = to_label, 
                    aes(label = to_print), 
                    show.legend = F) + 
    geom_hline(yintercept = PTHRESH, linetype = "dashed") +
    scale_alpha_manual(values = alpha_vals, guide = "none") + 
    scale_fill_manual(values = colpal_use) +
    scale_colour_manual(values = colpal_use) +
    scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "-log10(Enrichment PVAL)") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  return (res_plot)
}

cts_EUR_for_plot <- cts_EUR %>% 
  left_join(celltypes_tissue, by = c("Name" = "cts_name")) %>%
  mutate(neglogp = -log10(Coefficient_P_value),
         signif = ifelse(neglogp >= PTHRESH, "signif", "no")) %>%
  arrange(category, Name) 

# Order of names for plot
name_order <- unique(cts_EUR_for_plot$Name)
cts_EUR_for_plot$Name <- factor(cts_EUR_for_plot$Name,
                                levels = name_order)

lapply(c(HORMONES, INFERTILITY), function (hrinf) {
  ggsave(paste0(main_filepath, "/ovary_celltypes_heritability/plots/", 
                hrinf, "_EUR_cts.png"), 
         getEnrichmentPlot(cts_EUR_for_plot %>% filter(trait == hrinf)),
         height = 7, width = 7, units = "in")
})

ggsave(paste0(main_filepath, "/ovary_celltypes_heritability/plots/hormones_EUR_cts.png"), 
       getEnrichmentPlot(cts_EUR_for_plot %>% filter(trait %in% HORMONES)),
       height = 14, width = 14, units = "in")

ggsave(paste0(main_filepath, "/ovary_celltypes_heritability/plots/infertility_EUR_cts.png"), 
       getEnrichmentPlot(cts_EUR_for_plot %>% filter(trait %in% INFERTILITY)),
       height = 14, width = 14, units = "in")

# cts_all_anc_for_plot <- cts_all_anc %>%
#   left_join(celltypes_tissue, by = c("Name" = "cts_name")) %>%
#   mutate(neglogp = -log10(Coefficient_P_value),
#          signif = ifelse(neglogp >= PTHRESH, "signif", "no")) %>%
#   arrange(tissue_category, Name)
# 
# # Order of names for plot
# name_order <- unique(cts_all_anc_for_plot$Name)
# cts_all_anc_for_plot$Name <- factor(cts_all_anc_for_plot$Name,
#                                 levels = name_order)
# 
# ggsave(paste0(main_filepath, "/gtex_heritability/plots/hormones_all_anc_cts.png"),
#        getEnrichmentPlot(cts_all_anc_for_plot %>% filter(trait %in% HORMONES)),
#        height = 14, width = 14, units = "in")
# 
# ggsave(paste0(main_filepath, "/gtex_heritability/plots/infertility_all_anc_cts.png"),
#        getEnrichmentPlot(cts_all_anc_for_plot %>% filter(trait %in% INFERTILITY)),
#        height = 14, width = 14, units = "in")
# 

# For poster ----

PTHRESH <- 2.75 # This is -log10 for FDR < 5%
alpha_vals <- c(0.5, 1)
names(alpha_vals) <- c("no", "signif")

# Expand colour palette
colpal_use <- c("grey", "#D35C79")
names(colpal_use) <- c("no", "signif")

getEnrichmentBarPlot <- function (df) {
  res_plot <- ggplot(df, aes(x = Name, y = neglogp)) +
    geom_bar(aes(fill = signif), stat = "identity") +
    geom_hline(yintercept = PTHRESH, linetype = "dashed") +
    scale_fill_manual(values = colpal_use, guide = "none") +
    scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  return (res_plot)
}

cts_EUR_testo <- cts_EUR %>% filter(trait == "Testosterone_M") %>%
  left_join(celltypes_tissue, by = c("Name" = "cts_name")) %>%
  mutate(neglogp = -log10(Coefficient_P_value),
         signif = ifelse(neglogp >= PTHRESH, "signif", "no")) %>%
  arrange(tissue_category, desc(neglogp)) 

# Order of names for plot
name_order <- unique(cts_EUR_testo$Name)
cts_EUR_testo$Name <- factor(cts_EUR_testo$Name,
                                levels = name_order)

ggsave(paste0(main_filepath, "/gtex_heritability/plots/testo_M_EUR_cts_barplot.png"), 
       getEnrichmentBarPlot(cts_EUR_testo),
       height = 10, width = 25, units = "cm")

