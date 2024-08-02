# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/infertility_meta_mvp/debug"

# Read in data ----

dat <- read.table(paste0(mainpath, "/compare_mafs_laura_mvp_full.txt"), 
                  sep = "\t", header = T, stringsAsFactors = F)

dat <- dat %>%
  mutate(across(all_of(c("MAF_MVP", "MAF_old_meta", "FreqSE_new_meta")), as.numeric)) %>%
  filter(MAF_MVP > 0 & MAF_MVP < 1 & MAF_old_meta > 0 & MAF_old_meta < 1)

# MAF vs MAF plots coloured by FreqSE ----

# MAF bins by meta-analysis average MAF

log_bin_breaks <- ceiling(log10(min(dat$MAF_old_meta))):-1
bin_breaks <- c(0, 0.0001, 0.001, 0.01, 0.1, 1)
bin_labels <- c("<0.01%", "[0.01% - 0.1%)", "[0.1% - 1%)", "[1% - 10%)", ">=10%")
dat$MAF_bin <- cut(dat$MAF_old_meta,
                   breaks = bin_breaks, labels = bin_labels,
                   include.lowest = T)

mm_per_bin <- lapply(bin_labels, function (mb) {
  dat_for_plot <- dat %>%
    filter(MAF_bin == mb)
  
  mm_plot <- ggplot(dat_for_plot, 
                    aes(x = MAF_MVP, y = MAF_old_meta,
                        colour = FreqSE_new_meta)) +
    geom_point(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_gradient(low = "lightblue", high = "darkblue")
  
  return (mm_plot)
})

ggsave(paste0(mainpath, "/maf_vs_maf_plot.png"), 
       ggarrange(plotlist = mm_per_bin, nrow = 2, ncol = 3, 
                 common.legend = T))

# Correlation in MAF across bins

summ_dat <- dat %>%
  group_by(MAF_bin) %>%
  summarise(maf_correlation = cor(MAF_MVP, MAF_old_meta),
            avg_freqse = mean(FreqSE_new_meta))


