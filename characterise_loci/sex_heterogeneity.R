# Author: Samvida S. Venkatesh
# Date: 24/10/2022

library(argparse)
library(tidyverse)
theme_set(theme_bw())

# Get arguments ----

parser <- ArgumentParser()
parser$add_argument("--hormone", required = TRUE,
                    help = "Hormone to test sex-heterogeneity")
args <- parser$parse_args()

HORMONE <- args$hormone
logfile <- paste0("/well/lindgren/samvida/hormones_infertility/sex_heterogeneity/",
                  HORMONE, "_log.txt")

# Read data ----

filenames <- read.table("/well/lindgren/samvida/hormones_infertility/sex_heterogeneity/hormone_files_list.txt",
                        sep = "\t", header = F, stringsAsFactors = F)
colnames(filenames) <- c("strata", "fname", "N")

NUMERIC_COLS <- c("CHROM", "GENPOS", "Freq1", "A1FREQ", "BETA", "SE", "PVALUE")

dat_f <- read.table(filenames$fname[filenames$strata == paste0(HORMONE, "_F")], 
                    sep = "\t", header = T, stringsAsFactors = F)
dat_f <- dat_f %>%
  mutate(across(any_of(NUMERIC_COLS), as.numeric))
dat_f <- dat_f[, c(1:3, 5:7, 9:11)]
colnames(dat_f)[6:9] <- c("AF_Tested_F", "BETA_F", "SE_F", "PVALUE_F")

dat_m <- read.table(filenames$fname[filenames$strata == paste0(HORMONE, "_M")], 
                    sep = "\t", header = T, stringsAsFactors = F)
if (grepl("UKBB", filenames$fname[filenames$strata == paste0(HORMONE, "_M")])) {
  dat_m <- dat_m %>%
    mutate(across(any_of(NUMERIC_COLS), as.numeric))
  dat_m <- dat_m[, c(1:6, 8, 9, 12)]
  colnames(dat_m)[c(4:9)] <- c("Allele1", "Allele2",
                               "AF_Tested_M", "BETA_M", "SE_M", "PVALUE_M")
  
} else {
  dat_m <- dat_m %>%
    mutate(across(any_of(NUMERIC_COLS), as.numeric))
  dat_m <- dat_m[, c(1:3, 5:7, 9:11)]
  colnames(dat_m)[6:9] <- c("AF_Tested_M", "BETA_M", "SE_M", "PVALUE_M")
  
}

dat_sex_comb <- read.table(filenames$fname[filenames$strata == paste0(HORMONE, "_sex_comb")], 
                           sep = "\t", header = T, stringsAsFactors = F)
dat_sex_comb <- dat_sex_comb %>%
  mutate(across(any_of(NUMERIC_COLS), as.numeric))
dat_sex_comb <- dat_sex_comb[, c(1:3, 5:7, 9:11)]
colnames(dat_sex_comb)[6:9] <- c("AF_Tested_sex_comb", "BETA_sex_comb", 
                                 "SE_sex_comb", "PVALUE_sex_comb")

# Filter to only GWS SNPs ----

gws_f <- dat_f$ID[dat_f$PVALUE_F <= 5E-08]
gws_m <- dat_m$ID[dat_m$PVALUE_M <= 5E-08]
gws_sex_comb <- dat_sex_comb$ID[dat_sex_comb$PVALUE_sex_comb <= 5E-08]

retain_snps <- unique(c(gws_f, gws_m, gws_sex_comb))

dat_f <- dat_f %>% filter(ID %in% retain_snps)
dat_m <- dat_m %>% filter(ID %in% retain_snps)
dat_sex_comb <- dat_sex_comb %>% filter(ID %in% retain_snps)

# Combine data across sexes and check directional consistency ----

full_dat <- list(dat_f, dat_m, dat_sex_comb) %>% 
  reduce(inner_join)

full_dat <- full_dat %>% 
  mutate(sex_dirn_consistent = ifelse(BETA_F*BETA_M > 0, 
                                      "consistent", 
                                      "inconsistent"),
         het_zstat = (BETA_F - BETA_M)/sqrt(SE_F^2 + SE_M^2),
         het_pval = 2*pnorm(-abs(het_zstat), lower.tail = T))

write.table(full_dat,
            paste0("/well/lindgren/samvida/hormones_infertility/sex_heterogeneity/",
                   HORMONE, "_all_sexes_gws_results.txt"),
            sep = "\t", row.names = F, quote = F)

# Prune to get independent lead SNPs across analyses ----

PTHRESH <- 5E-08

all_snps <- full_dat %>%
  filter(PVALUE_F <= PTHRESH | PVALUE_M <= PTHRESH) %>%
  # filter(AF_Tested_F >= 0.01 & AF_Tested_F <= 0.99) %>%
  # filter(AF_Tested_M >= 0.01 & AF_Tested_M <= 0.99) %>%
  mutate(PVALUE = pmin(PVALUE_F, PVALUE_M))

all_snps <- split(all_snps, f = all_snps$CHROM)

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

# Bland-Altman plots ----

custom_three_diverge <- c("#D35C79","#009593", "#000000")
names(custom_three_diverge) <- c("female_only", "male_only", "both")

plotBlandAltman <- function (df) {
  # Decide colour
  for_plot <- df %>%
    mutate(beta_diff = BETA_F - BETA_M,
           beta_mean = (BETA_F + BETA_M)/2,
           point_type = ifelse(PVALUE_F <= 5E-8 & PVALUE_M <= 5E-8, "both",
                               ifelse(PVALUE_F <= 5E-8 & PVALUE_M > 5E-8, "female_only",
                                      ifelse(PVALUE_F > 5E-8 & PVALUE_M <= 5E-8, "male_only", NA)))) %>%
    filter(!is.na(point_type))
  
  # Get BA plot P-value for title
  ba_p_title <- t.test(for_plot$beta_diff, alternative = "two.sided",
                       mu = 0, paired = F, conf.level = 0.95)$p.val
  cat(paste0("Mean: ", signif(mean(for_plot$beta_diff), 3),
             ", P: ", signif(ba_p_title, 3)), "\n")
  
  res_plot <- ggplot(for_plot, aes(x = beta_mean, 
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

ggsave(paste0("/well/lindgren/samvida/hormones_infertility/sex_heterogeneity/",
              HORMONE, "_BA_plot.png"),
       plotBlandAltman(full_dat), 
       height = 5, width = 5, units = "cm")

# Scatter plots ----

plotScatter <- function (df) {
  # Decide colour
  for_plot <- df %>%
    mutate(LCI_M = BETA_M - 1.96*SE_M, UCI_M = BETA_M + 1.96*SE_M, 
           LCI_F = BETA_F - 1.96*SE_F, UCI_F = BETA_F + 1.96*SE_F, 
           point_type = ifelse(PVALUE_F <= 5E-8 & PVALUE_M <= 5E-8, "both",
                  ifelse(PVALUE_F <= 5E-8 & PVALUE_M > 5E-8, "female_only",
                         ifelse(PVALUE_F > 5E-8 & PVALUE_M <= 5E-8, "male_only", NA)))) %>%
    filter(!is.na(point_type))
  
  # Get R2 for title
  cat(cor(for_plot$BETA_M, for_plot$BETA_F))
  
  res_plot <- ggplot(for_plot, aes(x = BETA_M, y = BETA_F,
                                   fill = point_type, colour = point_type)) +
    geom_pointrange(aes(xmin = LCI_M, xmax = UCI_M),
                    size = 0.3, fatten = 2) +
    geom_pointrange(aes(ymin = LCI_F, ymax = UCI_F),
                    size = 0.3, fatten = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    #geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = custom_three_diverge, guide = "none") +
    scale_fill_manual(values = custom_three_diverge, guide = "none") +
    theme(axis.text = element_text(size = 6),
          axis.title = element_blank(),
          plot.title = element_text(size = 8))
  return (res_plot)
}

ggsave(paste0("/well/lindgren/samvida/hormones_infertility/sex_heterogeneity/",
              HORMONE, "_scatter_plot.png"),
       plotScatter(full_dat), 
       height = 5, width = 5, units = "cm")



