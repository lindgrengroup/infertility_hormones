# Author: Samvida S. Venkatesh
# Date: 26/09/2022

# R wrapper to submit pheno and strata-specific jobs
# Need to loop over all of the phenotypes and sex strata

main_filepath <- "/well/lindgren/samvida/hormones_infertility"
STANDARD_PHENO_NAMES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone")
STANDARD_SS_NAMES <- c("F", "M", "sex_comb")
STANDARD_ANC_NAMES <- c("all_anc", "EUR")

submission_script <- paste0(main_filepath, "/scripts/submit_filter_plot_meta_results.sh")

for (pheno in STANDARD_PHENO_NAMES) {
  for (anc in STANDARD_ANC_NAMES) {
    for (ss in STANDARD_SS_NAMES) {
      check_name <- paste0(main_filepath, "/meta_results_230613_no_ukb/", pheno,
                           "_", ss, "_", anc, "_1.tbl")
      # If file exists:
      if (file.exists(check_name)) {
        infile_name <- check_name
        logfile_name <- paste0(main_filepath, "/meta_results_230613_no_ukb/filtered/",
                               pheno, "_", ss, "_", anc, ".log")
        outfile_name <- paste0(main_filepath, "/meta_results_230613_no_ukb/filtered/", 
                               pheno, "_", ss, "_", anc, "_filtered.txt")
        outplot_dir <- paste0(main_filepath, "/meta_results_230613_no_ukb/plots/", 
                              pheno, "_", ss, "_", anc, "/")
        job_options <- paste0(
          "--export=",
          paste0(
            "inputFile=\"", infile_name, "\",",
            "logFile=\"", logfile_name, "\",",
            "outputFile=\"", outfile_name, "\",",
            "outPlotDir=\"", outplot_dir, "\""
          )
        )
        job_submission <- paste("sbatch", job_options, submission_script)
        system(job_submission)
        print(job_submission)
      }
    }
  }
}
