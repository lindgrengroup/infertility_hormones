# Author: Samvida S. Venkatesh
# Date: 24/10/2022

library(tidyverse)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility"

logfile <- paste0(mainpath, "/meta_results_230613_no_public/compare_public/log_file_240619.txt")

# Read data ----

NUMERIC_COLS <- c("Freq1", "Effect", "StdErr", "Pvalue")

STRATA <- read.table(paste0(mainpath, "/with_public_strata_list.txt"),
                     sep = "\t", header = F, stringsAsFactors = F)$V1

strata_dat <- lapply(STRATA, function (strata_name) {
  public_dat <- read.table(paste0(mainpath, "/meta_results_230613_no_public/compare_public/lead_snps_with_public_",
                                  strata_name, ".txt"),
                           sep = "\t", header = F, stringsAsFactors = F)
  colnames(public_dat) <- c("MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE",
                            "Effect", "StdErr", "Pvalue", "Direction", 
                            "HetISq", "HetChiSq", "HetDf", "HetPVal")
  public_dat <- public_dat %>%
    mutate(across(any_of(NUMERIC_COLS), as.numeric))
  public_dat <- public_dat[, c(1:4, 6:8)]
  colnames(public_dat)[4:7] <- c("AF_Tested_public", "BETA_public", "SE_public", "PVALUE_public")
  
  nopub_dat <- read.table(paste0(mainpath, "/meta_results_230613_no_public/compare_public/lead_snps_no_public_",
                                  strata_name, ".txt"),
                           sep = "\t", header = F, stringsAsFactors = F)
  colnames(nopub_dat) <- c("MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE",
                            "Effect", "StdErr", "Pvalue", "Direction", 
                            "HetISq", "HetChiSq", "HetDf", "HetPVal")
  nopub_dat <- nopub_dat %>%
    mutate(across(any_of(NUMERIC_COLS), as.numeric))
  nopub_dat <- nopub_dat[, c(1:4, 6:8)]
  colnames(nopub_dat)[4:7] <- c("AF_Tested_nopub", "BETA_nopub", "SE_nopub", "PVALUE_nopub")
  
  full_dat <- inner_join(public_dat, nopub_dat)
  return (full_dat)
})
names(strata_dat) <- STRATA

# Bland-Altman plots ----

plotBlandAltman <- function (df) {
  # Decide colour
  for_plot <- df %>%
    mutate(beta_diff = BETA_public - BETA_nopub,
           beta_mean = (BETA_public + BETA_nopub)/2) 
  
  # Get BA plot P-value for title
  ba_p_title <- t.test(for_plot$beta_diff, alternative = "two.sided",
                       mu = 0, paired = F, conf.level = 0.95)$p.val
  cat(paste0("Mean: ", signif(mean(for_plot$beta_diff), 3),
             ", P: ", signif(ba_p_title, 3)), "\n")
  
  res_plot <- ggplot(for_plot, aes(x = beta_mean, 
                                   y = beta_diff)) +
    geom_point() +
    #geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    #geom_abline(intercept = 0, slope = 1) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_blank(),
          plot.title = element_text(size = 8))
  return (res_plot)
}

# Report directional consistency ----

lapply(STRATA, function (strata_name) {
  summ_res <- strata_dat[[strata_name]] %>%
    mutate(signif_no_pub = PVALUE_nopub <= 5E-8,
           dirn_consistent = BETA_nopub*BETA_public > 0,
           het_zstat = (BETA_public - BETA_nopub)/sqrt(SE_public^2 + SE_nopub^2),
           het_pval = pnorm(het_zstat, 0, 1, lower.tail = T))
  
  sink(logfile, append = T)
  cat(paste0("STRATA: ", strata_name, "\n",
             "Number of GWS SNPs without public sumstats: ", sum(summ_res$signif_no_pub), "\n",
             "Number of directionally consistent SNPs: ", sum(summ_res$dirn_consistent), "\n",
             "Number of SNPs with heterogeneity PVAL < 0.05:", length(which(summ_res$het_pval < 0.05)), "\n"))
  sink()
  
  ggsave(paste0(mainpath, "/meta_results_230613_no_public/compare_public/plots/",
                strata_name, "_BA_plot.png"),
         plotBlandAltman(strata_dat[[strata_name]]), 
         height = 5, width = 5, units = "cm")
})




