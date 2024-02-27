# Author: Samvida S. Venkatesh
# Date: 17/04/2022

library(tidyverse)

mainpath <- "/well/lindgren/samvida/hormones_infertility/conditional_analysis"

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone") 
SEX_STRATA <- c("F", "M", "sex_comb")
ANCESTRY_GROUPS <- c("EUR", "all_anc")

sample_sizes <- read.table("/well/lindgren/samvida/hormones_infertility/meta_results_230613/sample_sizes.txt",
                           sep = "\t", header = T, stringsAsFactors = F)

SUBMISSION_SCRIPT <- "/well/lindgren/samvida/hormones_infertility/scripts/make_sumstats_for_gcta.sh"

# Apply across all strata ----

lapply(HORMONES, function (hr) {
  lapply(SEX_STRATA, function (sx) {
    lapply(ANCESTRY_GROUPS, function (anc_gp) {
      dir.create(paste0(mainpath, "/unreported_snps/",
                        hr, "_", sx, "_", anc_gp))
      job_options <- paste0(
        "--export=",
        paste0(
          "HORMONE=\"", hr, "\",",
          "SEX_STRATA=\"", sx, "\",",
          "ANC_GROUP=\"", anc_gp, "\",",
          "SAMPLE_SIZE=\"", sample_sizes$N[sample_sizes$strata == paste0(hr, "_", sx, "_", anc_gp)], "\""
        )
      )
      job_submission <- paste("sbatch", job_options, SUBMISSION_SCRIPT)
      system(job_submission)
      print(job_submission)
    })
  })
})
