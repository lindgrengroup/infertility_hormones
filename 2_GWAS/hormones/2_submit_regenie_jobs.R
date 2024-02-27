# Author: Samvida S. Venkatesh
# Date: 26/09/2022

# R wrapper to submit pheno and strata-specific jobs
# Need to loop over all of the phenotypes and sex strata

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone")
#ANC_GROUPS <- c("AFR", "EAS", "SAS")
ANC_GROUPS <- "non_finnish_EUR"
SEX_STRATA <- c("F", "M", "sex_comb")

#submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/scripts/regenie_step1.sh"
#submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/scripts/regenie_step2.sh"
submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/scripts/regenie_step2_xchr.sh"

for (hr in HORMONES) {
  for (anc in ANC_GROUPS) {
    for (ss in SEX_STRATA) {
      job_options <- paste0(
        "--export=",
        paste0(
          "HORMONE=\"", hr, "\",",
          "ANC_GROUP=\"", anc, "\",",
          "SEX_STRATA=\"", ss, "\""
        )
      )
      job_submission <- paste("sbatch", job_options, submission_script)
      system(job_submission)
      print(job_submission)
    }
  }
}

