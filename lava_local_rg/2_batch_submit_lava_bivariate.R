# Author: Samvida S. Venkatesh
# Date: 26/09/2022

# R wrapper to submit p1 x p2 jobs

main_filepath <- "/well/lindgren/samvida/hormones_infertility"

input_phenos <- read.table(paste0(main_filepath, "/lava_local/input_sumstats.txt"),
                           sep = "\t", header = T, stringsAsFactors = F)

# The first 13 phenotypes are the target (hormones, infertility, repro traits)
# everything else is obesity, which we don't want to run against other types of obesity

submission_script <- paste0(main_filepath, "/scripts/2_submit_lava_bivariate.sh")

for (i in 1:13) {
  for (j in (i+1):nrow(input_phenos)) {
    pheno1 <- input_phenos$phenotype[i]
    pheno2 <- input_phenos$phenotype[j]
    
    job_options <- paste0(
      "--export=",
      paste0(
        "pheno1=\"", pheno1, "\",",
        "pheno2=\"", pheno2, "\""
      )
    )
    job_submission <- paste("sbatch", job_options, submission_script)
    system(job_submission)
    print(job_submission)
    
  }
}

