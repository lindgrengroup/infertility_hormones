# Author: Samvida S. Venkatesh
# Date: 30/05/23

library(tidyverse)
library(TwoSampleMR)

# Read data ----

mainpath <- "/well/lindgren/samvida/hormones_infertility/two_sample_mr"

harmonised_dat <- readRDS(paste0(mainpath, "/harmonised_dat.rds"))

# Run Mendelian randomisation ----

# Remove Progesterone-F (0 SNPs)
KEEP_PAIRS <- paste0(rep(c("FSH_F", "LH_F", "Oestradiol_F", "Testosterone_F"),
                         each = 5), 
                     " x ", 
                     rep(c("female_infertility_analysis1", "female_infertility_analysis2", 
                           "female_infertility_analysis3", "female_infertility_analysis4", 
                           "female_infertility_analysis5"),
                         times = 4))
# Remove FSH-M (1 SNP), LH-M (1 SNP), Progesterone-M (0 SNPs) Oestradiol-M (1 SNP)
KEEP_PAIRS <- c(KEEP_PAIRS, "Testosterone_M x male_infertility")

mr_res <- mr(harmonised_dat)

# Test for pleiotropy - add MR Egger intercept
mr_pleio <- mr_pleiotropy_test(harmonised_dat)
colnames(mr_pleio) <- c("id.exposure", "id.outcome",
                        "outcome", "exposure",
                        "egger_intercept", "egger_intercept_se", "egger_intercept_pval")

sub_mr_res <- mr_res %>% 
  mutate(pair_tested = paste0(exposure, " x ", outcome)) %>%
  filter(pair_tested %in% KEEP_PAIRS)

# Add FDR adjustment within each method
sub_mr_res <- sub_mr_res %>%
  group_by(method) %>%
  mutate(pval.fdr = p.adjust(pval, method = "fdr"))

# Add pleiotropy test results
sub_mr_res <- inner_join(sub_mr_res,
                         mr_pleio)

write.table(sub_mr_res, paste0(mainpath, "/results_240621_with_egger_intercept.txt"),
            sep = "\t", row.names = F, quote = F)

