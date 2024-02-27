# Author: Samvida S. Venkatesh
# Date: 04/10/23

library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis3",
                   "female_infertility_analysis5")
REPRO_TRAITS <- c("dizygotic_twinning", 
                  "endometriosis", "heavy_menstrual_bleeding", 
                  "PCOS", "uterine_fibroids")
HORMONES <- c("TSH", "AMH", "FSH_F", "Testosterone_F")

# Read results ----

# reproductive trait results
rg_repro_results <- read.table("rg_phenome_wide/infertility_vs_repro_traits.txt",
                               sep = "\t", header = T, stringsAsFactors = F)

# within infertility strata
rg_infert_results <- read.table("rg_phenome_wide//infertility_vs_infertility.txt",
                                sep = "\t", header = T, stringsAsFactors = F)

# infertility vs hormones
rg_hr_results <- read.table("rg_hormones_infertility/EUR_genetic_correlations.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
rg_hr_results <- rg_hr_results %>%
  filter(trait_1 %in% HORMONES & trait_2 %in% INFERT_STRATA)

# phenome-wide
# Neale UKB manifest file with phenotype names
neale_ukb_manifest <- read.table("rg_phenome_wide/neale_phenotypes_manifest.txt",
                                 sep = "\t", header = T, stringsAsFactors = F,
                                 quote = "")
to_replace <- which(duplicated(neale_ukb_manifest$description))
neale_ukb_manifest$description[to_replace] <- paste0(neale_ukb_manifest$description[to_replace],
                                                     "_", neale_ukb_manifest$phenotype[to_replace])
sumstats_loc <- "/well/lindgren/resources/Neale_LDSC_all_phenos/"

rg_phenome <- lapply(INFERT_STRATA, function (infert) {
  res <- read.table(paste0("rg_phenome_wide/", infert, "_vs_all_phenos.r2"),
                    header = T, stringsAsFactors = F)
  res$trait_1 <- infert
  # Merge in p2 name from manifest
  res$p2 <- gsub(sumstats_loc, "", res$p2)
  res$trait_2 <- neale_ukb_manifest$description[match(res$p2,
                                                      neale_ukb_manifest$ldsc_sumstat_file)]
  return (res)
})
rg_phenome <- bind_rows(rg_phenome)

# Only keep phenos significant in at least one trait
# PTHRESH_FULL <- 0.05/(340*3)
# sig_rg_phenome <- rg_phenome %>% filter(p <= PTHRESH_FULL)
# phenos_keep <- unique(sig_rg_phenome$trait_2)
# rg_phenome <- rg_phenome %>% filter(trait_2 %in% phenos_keep)

# Keep chosen phenotypes
phenos_keep <- c("Ever had hysterectomy (womb removed)",
                 "Prospective memory result",
                 "Number of operations, self-reported",
                 "Tense / 'highly strung'",
                 "Age first had sexual intercourse",
                 "Fluid intelligence score",
                 "Body mass index (BMI)_23104_irnt",
                 "Waist circumference",
                 "Hip circumference",
                 "Body fat percentage",
                 "Comparative body size at age 10")
rg_phenome <- rg_phenome %>% filter(trait_2 %in% phenos_keep)

# Plot results ----

col_palette <- c("#D35C79", "#FFFFCC", "#009593")
names(col_palette) <- c("high", "mid", "low")
# overall min p-value for sizing points
max_psize <- ceiling(-log10(min(rg_repro_results$pval, rg_infert_results$pval,
                                rg_hr_results$pval, rg_phenome$p, na.rm = T)))

# Phenome-wide ----

# Order points 
to_order_rg_phenome <- rg_phenome %>%
  filter(trait_1 == "female_infertility_analysis1") %>%
  arrange(rg)
trait_2_levels <- to_order_rg_phenome$trait_2

pheno_wide_plot <- ggplot(rg_phenome, 
                          aes(x = factor(trait_1, levels = INFERT_STRATA), 
                              y = factor(trait_2, levels = trait_2_levels))) +        
  geom_tile(fill = "white", colour = "grey") +         
  geom_point(aes(fill = rg,
                 colour = rg, 
                 size = -log10(p))) +   
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79", limits = c(-1, 1), guide = "none") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79", limits = c(-1, 1), guide = "none") +     
  scale_size(range = c(1, 10), 
             limits = c(0, max_psize), guide = "none") +             
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank())

# ggsave("rg_phenome_wide/for_thesis_infertility_vs_all_phenos.png",
#        pheno_wide_plot,
#        units = "cm", height = 15, width = 5)

ggsave("rg_phenome_wide/for_manuscript_infertility_vs_all_phenos.png",
       pheno_wide_plot,
       units = "cm", height = length(phenos_keep), width = 5)

# Within infertility strata

infert_plot <- ggplot(rg_infert_results, 
                      aes(x = factor(trait_1, levels = INFERT_STRATA), 
                          y = factor(trait_2, levels = INFERT_STRATA))) +        
  geom_tile(fill = "white", colour = "grey") +         
  geom_point(aes(fill = rg,
                 colour = rg, 
                 size = -log10(pval))) +   
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79", limits = c(-1, 1), guide = "none") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79", limits = c(-1, 1), guide = "none") +     
  scale_size(range = c(1, 10), 
             limits = c(0, max_psize), guide = "none") +             
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank())
ggsave("rg_phenome_wide/for_thesis_infertility_vs_infertility.png",
       infert_plot,
       units = "cm", height = 3, width = 5)

# Hormones vs infertility
hr_plot <- ggplot(rg_hr_results, 
                  aes(x = factor(trait_2, levels = INFERT_STRATA), 
                      y = factor(trait_1, levels = HORMONES))) +        
  geom_tile(fill = "white", colour = "grey") +         
  geom_point(aes(fill = rg,
                 colour = rg, 
                 size = -log10(pval))) +   
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79", limits = c(-1, 1), guide = "none") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79", limits = c(-1, 1), guide = "none") +     
  scale_size(range = c(1, 10), 
             limits = c(0, max_psize), guide = "none") +             
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank())
ggsave("rg_phenome_wide/for_thesis_infertility_vs_hormones.png",
       hr_plot,
       units = "cm", height = 4, width = 5)


# Repro traits

repro_plot <- ggplot(rg_repro_results, 
                     aes(x = factor(trait_1, levels = INFERT_STRATA), 
                         y = factor(trait_2, levels = REPRO_TRAITS))) +        
  geom_tile(fill = "white", colour = "grey") +         
  geom_point(aes(fill = rg,
                 colour = rg, 
                 size = -log10(pval))) +   
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79", limits = c(-1, 1), guide = "none") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79", limits = c(-1, 1), guide = "none") +     
  scale_size(range = c(1, 10), 
             limits = c(0, max_psize), guide = "none") +             
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank())
ggsave("rg_phenome_wide/for_thesis_infertility_vs_repro.png",
       repro_plot,
       units = "cm", height = 5, width = 5)

