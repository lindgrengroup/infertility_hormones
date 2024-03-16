# Author: Samvida S. Venkatesh
# Date: 01/11/2021

library(tidyverse)
library(IRanges)

WINDOW_SIZE <- 10000 # 10kb around lead variant

# Read list of associations ----

# gwascat_assocns <- read.table("/well/lindgren/samvida/Resources/GWASCatalog/all_associations_231108.txt",
#                           sep = "\t", header = T, quote = "")
# 
# # Get chr and position
# gwascat_assocns_chrpos <- data.frame(chr = paste0("chr", gwascat_assocns$CHR_ID),
#                               pos = as.numeric(gwascat_assocns$CHR_POS))
# gwascat_assocns_chrpos <- gwascat_assocns_chrpos %>% distinct()
# 
# # Make sure there are no NAs
# gwascat_assocns_chrpos <- gwascat_assocns_chrpos[complete.cases(gwascat_assocns_chrpos), ]

gwascat_assocns <- read.table("/well/lindgren/samvida/Resources/GWASCatalog/gwascat_all_associations_hg38.bed",
                              sep = "\t", header = F, quote = "")
colnames(gwascat_assocns) <- c("chr", "pos")

# Create a list of non-overlapping ranges to see if variants fall within these ranges
intervals_check <- lapply(1:22, function (chrom) {
  sub_df <- gwascat_assocns %>%
    filter(chr == paste0("chr", chrom))
  intervals <- IRanges(start = sub_df$pos - WINDOW_SIZE,
                       end = sub_df$pos + WINDOW_SIZE)
  reduced_intervals <- reduce(intervals)
  return (reduced_intervals)
})

# Subset the scores to variants within +/- 10kb of GWASCat variants

sds_scores <- read.table("/well/lindgren/samvida/Resources/selection_scores_public/SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals_hg38.txt",
                         sep = " ", header = T, stringsAsFactors = F)
betascan_scores <- read.table("/well/lindgren/samvida/Resources/selection_scores_public/BetaScan2/hg38/all_B2std_hg38.txt",
                              sep = " ", header = T, stringsAsFactors = F)

# Go through selection scores and subset by position 
retain_sds_scores <- lapply(1:22, function (chrom) {
  sub_sds <- sds_scores %>% filter(CHR == paste0("chr", chrom))
  # Create an IRanges object for the positions
  number_ranges <- IRanges(start = sub_sds$POS, end = sub_sds$POS)
  # Find overlaps between number ranges and reduced intervals
  overlaps <- findOverlaps(number_ranges, intervals_check[[chrom]])
  # Extract the indices of numbers that fall within any interval
  pos_retain <- subjectHits(overlaps)
  return (sub_sds[unique(pos_retain), ])
})
retain_sds_scores <- bind_rows(retain_sds_scores)

# Go through selection scores and subset by position 
retain_beta_scores <- lapply(1:22, function (chrom) {
  sub_beta <- betascan_scores %>% filter(Chromosome == paste0("chr", chrom))
  # Create an IRanges object for the positions
  number_ranges <- IRanges(start = sub_beta$Position, end = sub_beta$Position)
  # Find overlaps between number ranges and reduced intervals
  overlaps <- findOverlaps(number_ranges, intervals_check[[chrom]])
  # Extract the indices of numbers that fall within any interval
  pos_retain <- subjectHits(overlaps)
  return (sub_beta[unique(pos_retain), ])
})
retain_beta_scores <- bind_rows(retain_beta_scores)

# Make sure integers as printed as integers and not in scientific
options(scipen = 999)
write.table(retain_sds_scores,
            "/well/lindgren/samvida/hormones_infertility/selection_analyses/data/gwascat_sds_scores.txt",
            sep = "\t", 
            row.names = F, quote = F, col.names = F)
options(scipen = 0)

options(scipen = 999)
write.table(retain_beta_scores,
            "/well/lindgren/samvida/hormones_infertility/selection_analyses/data/gwascat_betascan_scores.txt",
            sep = "\t", 
            row.names = F, quote = F, col.names = F)
options(scipen = 0)
