TRAITS <- c("FSH_F", "Testosterone_F", 
            "female_infertility_analysis1", "female_infertility_analysis3", "female_infertility_analysis5")
submission_script <- "/well/lindgren/samvida/hormones_infertility/scripts/calculate_rg.sh"

for (tt in TRAITS) {
  job_options <- paste0(
    "--export=",
    paste0(
      "STRATA=\"", tt, "\""
    )
  )
  job_submission <- paste("sbatch", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
