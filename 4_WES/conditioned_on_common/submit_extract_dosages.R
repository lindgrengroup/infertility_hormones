# Author: Samvida S. Venkatesh
# Date: 01/03/2022

# R wrapper to submit qctool jobs
# Loop over all lead variants in each parameter (across all strata)

library(tidyverse)

submission_script <- "/well/lindgren/samvida/hormones_infertility/scripts/3_get_dosages_common_variants.sh"

lead_variants <- read.table("/well/lindgren/samvida/hormones_infertility/exome_seq_results/conditioned_on_common/nearby_common_variants/unique_rsids_to_get.txt",
                            sep = "\t", header = F, stringsAsFactors = F)
colnames(lead_variants) <- c("SNP", "CHR", "POS")

auto_vars <- lead_variants %>% filter(CHR != 23)

for (vi in 1:nrow(auto_vars)) {
  snp <- auto_vars$SNP[vi]
  
  job_options <- paste0(
    "--export=",
    paste0(
      "CHR=\"", auto_vars$CHR[vi], "\",",
      "VARID=\"", snp, "\""
    )
  )
  job_submission <- paste("sbatch", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}

x_vars <- lead_variants %>% filter(CHR == 23)
submission_script <- "/well/lindgren/samvida/hormones_infertility/scripts/3_get_dosages_common_variants_xchr.sh"

for (vi in 1:nrow(x_vars)) {
  snp <- x_vars$SNP[vi]
  
  job_options <- paste0(
    "--export=",
    paste0(
      "VARID=\"", snp, "\""
    )
  )
  job_submission <- paste("sbatch", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
