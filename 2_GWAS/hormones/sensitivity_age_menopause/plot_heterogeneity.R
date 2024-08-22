# Author: Samvida S. Venkatesh
# Date: 24/10/2022

library(tidyverse)
theme_set(theme_bw())

# Read data ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/sensitivity_age_menopause"

HORMONES <- c("FSH", "LH")

NUMERIC_COLS <- c("GENPOS", "A1FREQ", 
                  "N_all", "BETA_all", "SE_all", "LOG10P_all",
                  "N_meno", "BETA_meno", "SE_meno", "LOG10P_meno")

# significant in UKB analyses (P<1E-5)
meno_dat <- lapply(HORMONES, function (hr) {
  df <- read.table(paste0(mainpath, "/signif_", hr, ".txt"),
                   sep = " ", header = T, stringsAsFactors = F)
  return (df)
})
names(meno_dat) <- HORMONES

# lead SNPs from meta-analysis
ma_lead_snps <- lapply(HORMONES, function (hr) {
  df <- read.table(paste0(mainpath, "/MA_lead_snps_", hr, ".txt"),
                   sep = " ", header = T, stringsAsFactors = F)
  return (df)
}) 
names(ma_lead_snps) <- HORMONES

# Wrangle data to check directional consistency ----

meno_dat <- lapply(HORMONES, function (hr) {
  df <- meno_dat[[hr]] %>%
    mutate(CHROM = ifelse(CHROM == 23, "X", CHROM),
           across(all_of(NUMERIC_COLS), as.numeric),
           PVALUE_all = 10^-LOG10P_all, PVALUE_meno = 10^-LOG10P_meno,
           dirn_consistent = ifelse(BETA_all*BETA_meno > 0, 
                                    "consistent", 
                                    "inconsistent"),
           het_zstat = (BETA_all - BETA_meno)/sqrt(SE_all^2 + SE_meno^2),
           het_pval = 2*pnorm(-abs(het_zstat), lower.tail = T))
  return (df)
})
names(meno_dat) <- HORMONES

ma_lead_snps <- lapply(HORMONES, function (hr) {
  df <- ma_lead_snps[[hr]] %>%
    mutate(CHROM = ifelse(CHROM == 23, "X", CHROM),
           across(all_of(NUMERIC_COLS), as.numeric),
           PVALUE_all = 10^-LOG10P_all, PVALUE_meno = 10^-LOG10P_meno,
           dirn_consistent = ifelse(BETA_all*BETA_meno > 0, 
                                    "consistent", 
                                    "inconsistent"),
           het_zstat = (BETA_all - BETA_meno)/sqrt(SE_all^2 + SE_meno^2),
           het_pval = 2*pnorm(-abs(het_zstat), lower.tail = T))
  
  write.table(df, paste0(mainpath, "/logs/MA_lead_snps_heterogeneity_",
                         hr, ".txt"), 
              sep = "\t", row.names = F, quote = F)
  return (df)
}) 
names(ma_lead_snps) <- HORMONES

# Prune to get independent lead SNPs across analyses ----

PTHRESH <- 1E-06
KB_PRUNE <- 500
WINDOW_SIZE <- KB_PRUNE*1000 # 500kb window

pruneSNPs <- function (df, dist_bp = 500000) {
  # Start with the first SNP in the dataframe
  old_top <- -Inf
  current_top <- df$GENPOS[1]
  # Check for SNP with lowest p-value in this region
  while (current_top != old_top) {
    old_top <- current_top
    current_dat <- df %>% filter(GENPOS >= (old_top - dist_bp) &
                                   GENPOS <= (old_top + dist_bp))
    current_top <- current_dat$GENPOS[which.min(current_dat$PVALUE)]
  }
  return (df[df$GENPOS == current_top, ])
}

meno_lead_snps <- lapply(HORMONES, function (hr) {
  all_snps <- meno_dat[[hr]] %>%
    mutate(PVALUE = pmin(PVALUE_all, PVALUE_meno))
  all_snps <- split(all_snps, f = all_snps$CHROM)
  
  # Apply to all genome-wide sig SNPs
  pruned_snps_full <- lapply(all_snps, function (cdf) {
    print(unique(cdf$CHROM))
    dat_to_prune <- cdf
    pruned_list <- c()
    dist_prune <- KB_PRUNE*1000
    # Repeat until there are no loci remaining
    while (nrow(dat_to_prune) > 0) {
      res_pruned <- pruneSNPs(dat_to_prune, dist_bp = WINDOW_SIZE)
      pruned_list <- c(pruned_list, res_pruned$ID)
      # Remove all SNPs within vicinity of pruned
      dat_to_prune <- dat_to_prune %>%
        filter(GENPOS <= (res_pruned$GENPOS - WINDOW_SIZE) |
                 GENPOS >= (res_pruned$GENPOS + WINDOW_SIZE))
    }
    # Results data
    pruned_dat <- cdf %>%
      filter(ID %in% pruned_list)
    return (pruned_dat)
  })
  pruned_snps_full <- bind_rows(pruned_snps_full)
  
  write.table(pruned_snps_full,
              paste0(mainpath, "/logs/ukb_independent_snps_heterogeneity_",
                     hr, ".txt"),
              sep = "\t", row.names = F, quote = F)
  
  return (pruned_snps_full)
})
names(meno_lead_snps) <- HORMONES

# Bland-Altman plots ----

custom_three_diverge <- c("#810098","#009139", "#000000")
names(custom_three_diverge) <- c("all_only", "meno_only", "both")

plotBlandAltman <- function (df) {
  # Decide colour
  for_plot <- df %>%
    mutate(beta_diff = BETA_all - BETA_meno,
           BETA_meno = (BETA_all + BETA_meno)/2,
           point_type = ifelse(PVALUE_all <= PTHRESH & PVALUE_meno <= PTHRESH, "both",
                               ifelse(PVALUE_all <= PTHRESH & PVALUE_meno > PTHRESH, "all_only",
                                      ifelse(PVALUE_all > PTHRESH & PVALUE_meno <= PTHRESH, "meno_only", NA)))) %>%
    filter(!is.na(point_type))
  
  # Get BA plot P-value for title
  ba_p_title <- t.test(for_plot$beta_diff, alternative = "two.sided",
                       mu = 0, paired = F, conf.level = 0.95)$p.val
  cat(paste0("Mean: ", signif(mean(for_plot$beta_diff), 3),
             ", P: ", signif(ba_p_title, 3)), "\n")
  
  res_plot <- ggplot(for_plot, aes(x = BETA_meno, 
                                   y = beta_diff,
                                   fill = point_type, colour = point_type)) +
    geom_point() +
    #geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    #geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = custom_three_diverge, guide = "none") +
    scale_fill_manual(values = custom_three_diverge, guide = "none") +
    theme(axis.text = element_text(size = 6),
          axis.title = element_blank(),
          plot.title = element_text(size = 8))
  return (res_plot)
}

# Scatter plots ----

plotScatter <- function (df) {
  # Decide colour
  for_plot <- df %>%
    mutate(LCI_meno = BETA_meno - 1.96*SE_meno, UCI_meno = BETA_meno + 1.96*SE_meno, 
           LCI_all = BETA_all - 1.96*SE_all, UCI_all = BETA_all + 1.96*SE_all, 
           point_type = ifelse(PVALUE_all <= PTHRESH & PVALUE_meno <= PTHRESH, "both",
                               ifelse(PVALUE_all <= PTHRESH & PVALUE_meno > PTHRESH, "all_only",
                                      ifelse(PVALUE_all > PTHRESH & PVALUE_meno <= PTHRESH, "meno_only", NA)))) %>%
    filter(!is.na(point_type))
  
  # Get R2 for title
  cat(paste0("R2: ", cor(for_plot$BETA_meno, for_plot$BETA_all)), "\n")
  
  # Get axis limits (same for x and y)
  min_axis <- min(c(for_plot$LCI_meno, for_plot$LCI_all))
  max_axis <- max(c(for_plot$UCI_meno, for_plot$UCI_all))
  
  res_plot <- ggplot(for_plot, aes(x = BETA_meno, y = BETA_all,
                                   fill = point_type, colour = point_type)) +
    geom_pointrange(aes(xmin = LCI_meno, xmax = UCI_meno),
                    size = 0.3, fatten = 2) +
    geom_pointrange(aes(ymin = LCI_all, ymax = UCI_all),
                    size = 0.3, fatten = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    #geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = custom_three_diverge, guide = "none") +
    scale_fill_manual(values = custom_three_diverge, guide = "none") +
    scale_x_continuous(limits = c(min_axis, max_axis)) +
    scale_y_continuous(limits = c(min_axis, max_axis)) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_blank(),
          plot.title = element_text(size = 8))
  return (res_plot)
}

lapply(HORMONES, function (hr) {
  ggsave(paste0(mainpath, "/plots/", hr, "_ukb_snps_BA_plot.png"),
         plotBlandAltman(meno_lead_snps[[hr]]), 
         height = 5, width = 5, units = "cm")
  
  ggsave(paste0(mainpath, "/plots/", hr, "_ukb_snps_scatter_plot.png"),
         plotScatter(meno_lead_snps[[hr]]), 
         height = 5, width = 5, units = "cm")
  
  ggsave(paste0(mainpath, "/plots/", hr, "_MA_lead_snps_BA_plot.png"),
         plotBlandAltman(ma_lead_snps[[hr]]), 
         height = 5, width = 5, units = "cm")
  
  ggsave(paste0(mainpath, "/plots/", hr, "_MA_lead_snps_scatter_plot.png"),
         plotScatter(ma_lead_snps[[hr]]), 
         height = 5, width = 5, units = "cm")
})


# Mean: 0.0225, P: 0.779
# R2: 0.793669882917108
# Mean: -0.0382, P: 0.249
# R2: 0.992668897499758
# Mean: -0.134, P: 0.427
# R2: 0.805529556479772
# Mean: 0.0289, P: 0.606
# R2: 1



