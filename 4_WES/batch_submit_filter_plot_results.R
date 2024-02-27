# Author: Samvida S. Venkatesh
# Date: 26/09/2022

# R wrapper to submit pheno and strata-specific jobs
# Need to loop over all of the phenotypes and sex strata

STRATA <- c("FSH_F", "FSH_M", "FSH_sex_comb",
            "LH_F", "LH_M", "LH_sex_comb",
            "Oestradiol_F", "Oestradiol_M", "Oestradiol_sex_comb",
            "Testosterone_F", "Testosterone_M", "Testosterone_sex_comb",
            "female_infertility_binary", "idiop_infertility_exclusion_binary",
            "male_infertility_binary")

STRATA <- c("female_infertility_binary", "idiop_infertility_exclusion_binary",
            "male_infertility_binary")

submission_script <- "/well/lindgren/samvida/hormones_infertility/exome_seq_results/scripts/2_submit_filter_plot_results.sh"

for (str in STRATA) {
  job_options <- paste0(
    "--export=",
    paste0(
      "STRATA=\"", str, "\""
    )
  )
  job_submission <- paste("sbatch", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
