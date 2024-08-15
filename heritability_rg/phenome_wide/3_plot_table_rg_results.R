# Author: Samvida S. Venkatesh
# Date: 04/10/23

library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/rg_phenome_wide"
HORMONE_STRATA <- c("FSH_F", "Testosterone_F", "Testosterone_M")
INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis3",
                   "female_infertility_analysis5")

# Fairly conservative P-threshold adjusting for 500 tests 
PTHRESH <- 0.05/500

# Read results ----

# Neale UKB manifest file with phenotype names
neale_ukb_manifest <- read.table(paste0(mainpath, "/neale_ukb_pheno_manifest.txt"),
                                 sep = "\t", header = T, stringsAsFactors = F,
                                 quote = "")
to_replace <- which(duplicated(neale_ukb_manifest$description))
neale_ukb_manifest$description[to_replace] <- paste0(neale_ukb_manifest$description[to_replace],
                                                     "_", neale_ukb_manifest$phenotype[to_replace])
sumstats_loc <- "/well/lindgren/resources/Neale_LDSC_all_phenos/"

# rG results
rg_results <- lapply(c(HORMONE_STRATA, INFERT_STRATA), function (tt) {
  res <- read.table(paste0(mainpath, "/", tt, "_vs_all_phenos.r2"),
                    header = T, stringsAsFactors = F)
  res$p1 <- tt
  # Merge in p2 name from manifest
  res$p2 <- gsub(sumstats_loc, "", res$p2)
  res$p2_desc <- neale_ukb_manifest$description[match(res$p2,
                                                     neale_ukb_manifest$ldsc_sumstat_file)]
  return (res)
})
rg_results <- bind_rows(rg_results)

write_sig <- rg_results %>%
  #filter(p <= PTHRESH) %>%
  arrange(p1, p)

write.table(write_sig, paste0(mainpath, "/all_results_all_phenos.txt"),
             sep = "\t", quote = F, row.names = F)

# Plot significant results ----

# Hormones plot
sig_hr_phenos <- rg_results %>% filter(p1 %in% HORMONE_STRATA &
                                         p <= PTHRESH)
sig_hr_phenos <- unique(sig_hr_phenos$p2)

plot_hr <- rg_results %>%
  filter(p1 %in% HORMONE_STRATA & p2 %in% sig_hr_phenos) %>%
  mutate(p2_for_plot = stringr::str_trunc(p2_desc, 30)) %>%
  arrange(rg)

pheno_order <- unique(plot_hr$p2_for_plot)

col_palette <- c("#D35C79", "#FFFFCC", "#009593")
names(col_palette) <- c("high", "mid", "low")

hr_plot <- ggplot(plot_hr, aes(x = factor(p1, levels = HORMONE_STRATA), 
                    y = factor(p2_for_plot, levels = pheno_order))) +        
  geom_tile(fill = "white", colour = "grey") +         
  geom_point(aes(fill = rg,
                 colour = rg, 
                 size = -log10(p))) +   
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79") +     
  scale_size(range = c(1, 5)) +             
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0(mainpath, "/hormones_rg_plot.pdf"),
       hr_plot)

# Infertility plot
sig_infert_phenos <- rg_results %>% filter(p1 %in% INFERT_STRATA &
                                         p <= PTHRESH)
sig_infert_phenos <- unique(sig_infert_phenos$p2)

plot_infert <- rg_results %>%
  filter(p1 %in% INFERT_STRATA & p2 %in% sig_infert_phenos) %>%
  mutate(p2_for_plot = stringr::str_trunc(p2_desc, 30)) %>%
  arrange(rg)

pheno_order <- unique(plot_infert$p2_for_plot)

col_palette <- c("#D35C79", "#FFFFCC", "#009593")
names(col_palette) <- c("high", "mid", "low")

infert_plot <- ggplot(plot_infert, aes(x = factor(p1, levels = INFERT_STRATA), 
                               y = factor(p2_for_plot, levels = pheno_order))) +        
  geom_tile(fill = "white", colour = "grey") +         
  geom_point(aes(fill = rg,
                 colour = rg, 
                 size = -log10(p))) +   
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79") +
  scale_size(range = c(1, 5)) +             
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0(mainpath, "/infertility_rg_plot.pdf"),
       infert_plot)

## For manuscript main figure, pick certain phenotypes

PHENOS_PLOT <- c("Ever had hysterectomy (womb removed)",
                 "Prospective memory result",
                 "Number of operations, self-reported",
                 "Tense / 'highly strung'",
                 "Age first had sexual intercourse",
                 "Fluid intelligence score",
                 "Body mass index (BMI)_23104_irnt",
                 "Waist circumference",
                 "Hip circumference",
                 "Body fat percentage")

sig_infert_phenos <- rg_results %>% filter(p1 %in% INFERT_STRATA &
                                             p2_desc %in% PHENOS_PLOT)
plot_infert <- sig_infert_phenos %>%
  mutate(p2_for_plot = stringr::str_trunc(p2_desc, 30)) %>%
  arrange(rg)

pheno_order <- unique(plot_infert$p2_for_plot)

col_palette <- c("#D35C79", "#FFFFCC", "#009593")
names(col_palette) <- c("high", "mid", "low")

infert_plot <- ggplot(plot_infert, aes(x = factor(p1, levels = INFERT_STRATA), 
                                       y = factor(p2_for_plot, levels = pheno_order))) +        
  geom_tile(fill = "white", colour = "grey") +         
  geom_point(aes(fill = rg,
                 colour = rg, 
                 size = -log10(p))) +   
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79") +
  scale_size(range = c(1, 5)) +             
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0(mainpath, "/infertility_rg_plot.pdf"),
       infert_plot)


