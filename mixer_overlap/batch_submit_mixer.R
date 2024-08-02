# Author: Samvida S. Venkatesh
# Date: 26/09/2022

# R wrapper to submit p1 x p2 jobs

main_filepath <- "/well/lindgren/samvida/hormones_infertility"

input_phenos <- read.table(paste0(main_filepath, "/mixer_overlap/input_sumstats.txt"),
                           sep = "\t", header = T, stringsAsFactors = F)

# Run fit 1 
submission_script <- paste0(main_filepath, "/scripts/run_mixer_fit1.sh")

for (i in 1:nrow(input_phenos)) {
    pheno <- input_phenos$phenotype[i]
    pheno_loc <- input_phenos$filename[i]
    
    job_options <- paste0(
      "--export=",
      paste0(
        "pheno=\"", pheno, "\",",
        "pheno_loc=\"", pheno_loc, "\""
      )
    )
    job_submission <- paste("sbatch", job_options, submission_script)
    system(job_submission)
    print(job_submission)
    
}

# Run fit 2

# The first 3 phenotypes (infertility) are the target

submission_script <- paste0(main_filepath, "/scripts/run_mixer_fit2.sh")

for (i in 1:3) {
  for (j in 4:nrow(input_phenos)) {
    pheno1 <- input_phenos$phenotype[i]
    pheno2 <- input_phenos$phenotype[j]
    
    pheno1_loc <- input_phenos$filename[i]
    pheno2_loc <- input_phenos$filename[j]
    
    job_options <- paste0(
      "--export=",
      paste0(
        "pheno1=\"", pheno1, "\",",
        "pheno2=\"", pheno2, "\",",
        "pheno1_loc=\"", pheno1_loc, "\",",
        "pheno2_loc=\"", pheno2_loc, "\""
      )
    )
    job_submission <- paste("sbatch", job_options, submission_script)
    system(job_submission)
    print(job_submission)
    
  }
}
