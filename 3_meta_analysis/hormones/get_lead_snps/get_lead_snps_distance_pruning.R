# Author: Samvida S. Venkatesh
# Date: 17/04/2022

library(argparse)
library(tidyverse)

mainpath <- "/well/lindgren/samvida/hormones_infertility/meta_results_230613/lead_snps"

parser <- ArgumentParser()
parser$add_argument("--hormone", required=TRUE)
parser$add_argument("--sex_strata", required=TRUE)
parser$add_argument("--anc_group", required=TRUE)
args <- parser$parse_args()

HR <- args$hormone
SS <- args$sex_strata
ANC <- args$anc_group

# GWS SNPs ----

all_gws_snps <- read.table(paste0(mainpath, "/gws_variants_",
                                  HR, "_", SS, "_", ANC, "_with_rsids.txt"),
                           header = T, stringsAsFactors = F)
varids_1kg <- read.table(paste0(mainpath, "/plink_log/retain_1000G_variants.txt"),
                           header = F, stringsAsFactors = F)$V1

tmp <- all_gws_snps %>%
  mutate(poss_id_1kg = paste0(CHROM, ":", GENPOS, ":", Allele1, ":", Allele2),
         keep_no_flip = poss_id_1kg %in% varids_1kg,
         flipped_id_1kg = paste0(CHROM, ":", GENPOS, ":", Allele2, ":", Allele1),
         keep_flip = flipped_id_1kg %in% varids_1kg) 

# Create the 1kg varid for classification downstream
gws_snps_1kg <- all_gws_snps[tmp$keep_no_flip, ]
gws_snps_1kg <- gws_snps_1kg %>%
  mutate(VARID_1KG = paste0(CHROM, ":", GENPOS, ":", Allele1, ":", Allele2))
to_bind <- all_gws_snps[tmp$keep_flip, ]
to_bind <- to_bind %>%
  mutate(VARID_1KG = paste0(CHROM, ":", GENPOS, ":", Allele2, ":", Allele1))

gws_snps_1kg <- bind_rows(gws_snps_1kg, to_bind)

# Split by chromosome for pruning (can combine again later)
all_gws_snps <- split(all_gws_snps, f = all_gws_snps$CHROM)
gws_snps_1kg <- split(gws_snps_1kg, f = gws_snps_1kg$CHROM)

# Distance-based pruning of SNPs ----

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
pruned_snps_full <- lapply(all_gws_snps, function (cdf) {
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
  pruned_dat$hormone <- HR
  pruned_dat$sex_strata <- SS
  pruned_dat$ancestry <- ANC
  return (pruned_dat)
})
pruned_snps_full <- bind_rows(pruned_snps_full)

# Write to table
write.table(pruned_snps_full,
            paste0(mainpath, "/all_lead_snp_sumstats_",
                   HR, "_", SS, "_", ANC, "_with_rsids.txt"),
            row.names = F, quote = F, sep = "\t")

# Apply to only 1KG SNPs

pruned_snps_1kg <- lapply(gws_snps_1kg, function (cdf) {
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
  pruned_dat$hormone <- HR
  pruned_dat$sex_strata <- SS
  pruned_dat$ancestry <- ANC
  return (pruned_dat)
})
pruned_snps_1kg <- bind_rows(pruned_snps_1kg)

# Write to table
write.table(pruned_snps_1kg,
            paste0(mainpath, "/1kg_lead_snp_sumstats_",
                   HR, "_", SS, "_", ANC, "_with_rsids.txt"),
            row.names = F, quote = F, sep = "\t")

