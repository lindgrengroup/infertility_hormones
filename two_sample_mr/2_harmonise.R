# Author: Samvida S. Venkatesh
# Date: 30/05/23

library(tidyverse)
library(TwoSampleMR)

# Read data ----

mainpath <- "/well/lindgren/samvida/hormones_infertility/two_sample_mr"

HORMONES <- c("FSH_F", "LH_F", "Oestradiol_F", "Progesterone_F",
              "Testosterone_F", "Testosterone_M")
INFERTILITY <- c("female_infertility_analysis1",
                 "female_infertility_analysis2",
                 "female_infertility_analysis3",
                 "female_infertility_analysis4",
                 "female_infertility_analysis5",
                 "male_infertility")

# Ensure instruments exist
instrument_filenames <- paste0(mainpath, "/exposures/all_lead_snp_sumstats_", HORMONES,
                               "_EUR_with_rsids.txt")
names(instrument_filenames) <- HORMONES
HORMONES <- HORMONES[file.exists(instrument_filenames)]
# Ensure outcomes exis
outcome_filenames <- paste0(mainpath, "/outcomes/instrument_outcomes_for_", 
                            INFERTILITY, "_EUR.txt")
names(outcome_filenames) <- INFERTILITY
INFERTILITY <- INFERTILITY[file.exists(outcome_filenames)]

# Read exposures
instruments <- lapply(HORMONES, function (hr) {
  res <- read_exposure_data(filename = instrument_filenames[hr],
                            sep = "\t",
                            snp_col = "RSID",
                            beta_col = "BETA",
                            se_col = "SE",
                            effect_allele_col = "Allele1",
                            other_allele_col = "Allele2",
                            eaf_col = "Freq1",
                            pval_col = "PVALUE")
  res$exposure <- hr
  return (res)
})
instruments <- bind_rows(instruments)

# Read outcomes

# When doing this the first time, may have to manually check for duplicates
# and make sure alleles match
# rs564372300 T/G

outcomes <- lapply(INFERTILITY, function (i) {
  res <- read_outcome_data(filename = outcome_filenames[i],
                           sep = "\t",
                           snp_col = "RSID",
                           beta_col = "BETA",
                           se_col = "SE",
                           effect_allele_col = "ALLELE1",
                           other_allele_col = "ALLELE2",
                           eaf_col = "FREQ1",
                           pval_col = "PVALUE")
  res$outcome <- i
  return (res)
})
outcomes <- bind_rows(outcomes)

# Harmonise data ----

# We don't need to infer alleles because we've performed the GWASs and
# we know the strand (i.e. don't try to flip strands for ambiguous variants)
harmonised_dat <- harmonise_data(instruments, outcomes, action = 1)

saveRDS(harmonised_dat, paste0(mainpath, "/harmonised_dat.rds"))


