# Author: Samvida S. Venkatesh
# Date: 02/10/23

library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/selection_analyses"

INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis2",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5",
                   "male_infertility")

# Read lead SNPs ----

lead_snps <- read.table(paste0(mainpath, "/all_hormones_infertility_lead_variants.txt"),
                        sep = "\t", header = T, stringsAsFactors = F)

# Read selection scores ----

sds_scores <- read.table("/well/lindgren/samvida/Resources/selection_scores_public/SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals_hg38.txt",
                         sep = " ", header = T, stringsAsFactors = F)
SDS_LOW <- quantile(sds_scores$SDS, 0.025)
SDS_HIGH <- quantile(sds_scores$SDS, 0.975)

betascan_scores <- read.table("/well/lindgren/samvida/Resources/selection_scores_public/BetaScan2/hg38/all_B2std_hg38.txt",
                              sep = " ", header = T, stringsAsFactors = F)
BETA_HIGH <- quantile(betascan_scores$Beta2_std, 0.95)

# Read scores for GWASCatalog regions only ----

sds_gwascat <- read.table("/well/lindgren/samvida/hormones_infertility/selection_analyses/data/gwascat_sds_scores.txt",
                          sep = "\t", header = F, stringsAsFactors = F)
colnames(sds_gwascat) <- colnames(sds_scores)
sds_perc <- ecdf(sort(sds_gwascat$SDS))

betascan_gwascat <- read.table("/well/lindgren/samvida/hormones_infertility/selection_analyses/data/gwascat_betascan_scores.txt",
                               sep = "\t", header = F, stringsAsFactors = F)
colnames(betascan_gwascat) <- colnames(betascan_scores)
betascan_perc <- ecdf(sort(betascan_gwascat$Beta2_std))

# Function to extract scores for all variants within 10 kb of lead variant ----

WINDOW_SIZE <- 10000 # bp
getWindowScores <- function (chr_lead, pos_lead) {
  sds_window <- sds_scores %>% filter(CHR == paste0("chr", chr_lead) &
                                        POS >= pos_lead - WINDOW_SIZE &
                                        POS <= pos_lead + WINDOW_SIZE) %>%
    mutate(status = ifelse(SDS >= SDS_HIGH | SDS <= SDS_LOW, "signif", "non_sig")) %>%
    rename(BP_pos = POS, score = SDS)
  
  beta_window <- betascan_scores %>% filter(Chromosome == paste0("chr", chr_lead) &
                                              Position >= pos_lead - WINDOW_SIZE &
                                              Position <= pos_lead + WINDOW_SIZE) %>%
    mutate(status = ifelse(Beta2_std >= BETA_HIGH, "signif", "non_sig")) %>%
    rename(BP_pos = Position, score = Beta2_std)
  
  return (list(sds = sds_window,
               betascan = beta_window))
}

# Apply across all infertility traits ----

infert_selection_scores <- lapply(INFERT_STRATA, function (infert) {
  snps_to_plot <- lead_snps %>% filter(trait == infert)
  cat(paste0("Trait: ", infert), "\n")
  if (nrow(snps_to_plot) > 0) {
    selection_scores <- lapply(1:nrow(snps_to_plot), function (i) {
      return (getWindowScores(chr_lead = snps_to_plot$chrom[i],
                              pos_lead = snps_to_plot$pos_hg38[i]))
    })
    names(selection_scores) <- snps_to_plot$rsid
    
    # Get summary of min/max SDS and max Beta scores
    summ_dat <- lapply(snps_to_plot$rsid, function (i) { 
      sds_df <- selection_scores[[i]]$sds
      beta_df <- selection_scores[[i]]$betascan
      return (data.frame(rsid = i,
                         n_SDS = nrow(sds_df),
                         min_SDS = min(sds_df$score),
                         min_SDS_perc = sds_perc(min(sds_df$score)),
                         max_SDS = max(sds_df$score),
                         max_SDS_perc = sds_perc(max(sds_df$score)),
                         n_beta = nrow(beta_df),
                         max_beta = max(beta_df$score),
                         max_beta_perc = betascan_perc(max(beta_df$score))))
    })
    summ_dat <- bind_rows(summ_dat)
    summ_dat$trait <- infert
    } else summ_dat <- NULL
  return (summ_dat)
})
infert_selection_scores <- bind_rows(infert_selection_scores)
write.table(infert_selection_scores,
            paste0(mainpath, "/results/infertility_loci_selection_scores.txt"),
            sep = "\t", row.names = F, quote = F)

# Format for thesis table
infert_tt <- infert_selection_scores %>%
  mutate(write_min_sds = paste0(signif(min_SDS, 3), "(", signif(min_SDS_perc*100, 3), "%)"),
         write_max_sds = paste0(signif(max_SDS, 3), "(", signif(max_SDS_perc*100, 3), "%)"),
         write_max_beta = paste0(signif(max_beta, 3), "(", signif(max_beta_perc*100, 3), "%)")) %>%
  select(all_of(c("trait", "rsid",
                  "n_SDS", "write_min_sds", "write_max_sds",
                  "n_beta", "write_max_beta")))
write.table(infert_tt,
            paste0(mainpath, "/results/infertility_loci_selection_scores_for_tt.txt"),
            sep = "\t", row.names = F, quote = F)
