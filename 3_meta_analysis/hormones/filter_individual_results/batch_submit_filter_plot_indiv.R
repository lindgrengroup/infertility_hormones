# Author: Samvida S. Venkatesh
# Date: 26/09/2022

# R wrapper to submit pheno and strata-specific jobs
# Need to loop over all of the phenotypes and sex strata

main_filepath <- "/well/lindgren/samvida/hormones_infertility"
STANDARD_PHENO_NAMES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone")
STANDARD_ANC_NAMES <- c("AFR", "EAS", "SAS", "non_finnish_EUR")
STANDARD_SS_NAMES <- c("F", "M", "sex_comb")

# ALSPAC ----

submission_script <- paste0(main_filepath, "/scripts/submit_filter_plot_indiv_BOLT.sh")

ALSPAC_phenos <- read.table(paste0(main_filepath, "/data/ALSPAC/pheno_list.txt"),
                            sep = "\t", header = F, stringsAsFactors = F)$V1

# Remap names so that they are standardised across all analyses
names(ALSPAC_phenos) <- c("FSH", "LH", "Testosterone")

for (pheno in c("FSH", "LH", "Testosterone")) {
  infile_name <- paste0(main_filepath, "/data/ALSPAC/lifted/", 
                        ALSPAC_phenos[pheno], "_hg38.txt")
  # Not all pheno-sexstrata combinations have input, so check that the file exists
  if (file.exists(infile_name)) {
    logfile_name <- paste0(main_filepath, "/data/ALSPAC/filtered/", pheno, "_F.log")
    outfile_name <- paste0(main_filepath, "/data/ALSPAC/filtered/", pheno, "_F_filtered.txt")
    outplot_dir <- paste0(main_filepath, "/plots/ALSPAC/", pheno, "_F/")
    
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

# Estonian Biobank ----

submission_script <- paste0(main_filepath, "/scripts/submit_filter_plot_indiv_regenie.sh")

EstBB_phenos <- read.table(paste0(main_filepath, "/data/EstBB/pheno_list.txt"),
                           sep = "\t", header = F, stringsAsFactors = F)$V1
EstBB_sexstrata <- c("female", "male", "combined")

# Remap names so that they are standardised across all analyses
names(EstBB_phenos) <- STANDARD_PHENO_NAMES
names(EstBB_sexstrata) <- STANDARD_SS_NAMES

for (pheno in STANDARD_PHENO_NAMES) {
  for (ss in STANDARD_SS_NAMES) {
    infile_name <- paste0(main_filepath, "/data/EstBB/lifted/", 
                          EstBB_phenos[pheno], "_", EstBB_sexstrata[ss], 
                          "_hg38.txt")
    # Not all pheno-sexstrata combinations have input, so check that the file exists
    if (file.exists(infile_name)) {
      logfile_name <- paste0(main_filepath, "/data/EstBB/filtered/", pheno, "_", ss, 
                             ".log")
      outfile_name <- paste0(main_filepath, "/data/EstBB/filtered/", pheno, "_", ss, 
                             "_filtered.txt")
      outplot_dir <- paste0(main_filepath, "/plots/EstBB/", pheno, "_", ss, "/")
      
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

# Genes & Health ----

submission_script <- paste0(main_filepath, "/scripts/submit_filter_plot_indiv_regenie.sh")

GH_phenos <- read.table(paste0(main_filepath, "/data/GenesHealth/pheno_list.txt"),
                        sep = "\t", header = F, stringsAsFactors = F)$V1
GH_sexstrata <- c("females", "males", "all")

# Remap names so that they are standardised across all analyses
names(GH_phenos) <- STANDARD_PHENO_NAMES
names(GH_sexstrata) <- STANDARD_SS_NAMES

for (pheno in STANDARD_PHENO_NAMES) {
  for (ss in STANDARD_SS_NAMES) {
    infile_name <- paste0(main_filepath, "/data/GenesHealth/genes_and_health_quant_traits_",
                          GH_sexstrata[ss], "_SAS_Jacobs_10_11_2022_", 
                          GH_phenos[pheno], ".regenie")
    # Not all pheno-sexstrata combinations have input, so check that the file exists
    if (file.exists(infile_name)) {
      logfile_name <- paste0(main_filepath, "/data/GenesHealth/filtered/", pheno, "_", ss, 
                             ".log")
      outfile_name <- paste0(main_filepath, "/data/GenesHealth/filtered/", pheno, "_", ss, 
                             "_filtered.txt")
      outplot_dir <- paste0(main_filepath, "/plots/GenesHealth/", pheno, "_", ss, "/")
      
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

# UK Biobank ----

submission_script <- paste0(main_filepath, "/scripts/submit_filter_plot_indiv_regenie.sh")
#submission_script <- paste0(main_filepath, "/scripts/tmp.sh")

for (pheno in STANDARD_PHENO_NAMES) {
  for (anc in STANDARD_ANC_NAMES) {
    for (ss in STANDARD_SS_NAMES) {
      infile_name <- paste0(main_filepath, "/data/UKBB/lifted/", 
                            pheno, "_", anc, "_", ss, "_hg38.txt")
      # Not all pheno-anc-sexstrata combinations have input, so check that the file exists
      if (file.exists(infile_name)) {
        logfile_name <- paste0(main_filepath, "/data/UKBB/filtered/", pheno, "_", anc, "_", ss, 
                               ".log")
        outfile_name <- paste0(main_filepath, "/data/UKBB/filtered/", pheno, "_", anc, "_", ss, 
                               "_filtered.txt")
        outplot_dir <- paste0(main_filepath, "/plots/UKBB/", pheno, "_", anc, "_", ss, "/")
        
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

# deCODE ----

submission_script <- paste0(main_filepath, "/scripts/submit_filter_plot_deCode.sh")

deCODE_phenos <- read.table(paste0(main_filepath, "/data/deCode/pheno_list.txt"),
                        sep = "\t", header = F, stringsAsFactors = F)$V1
deCODE_sexstrata <- c("female", "male", "combined")

# Remap names so that they are standardised across all analyses
names(deCODE_phenos) <- STANDARD_PHENO_NAMES
names(deCODE_sexstrata) <- STANDARD_SS_NAMES

for (pheno in STANDARD_PHENO_NAMES) {
  for (ss in STANDARD_SS_NAMES) {
    infile_name <- paste0(main_filepath, "/data/deCode/",
                          deCODE_phenos[pheno], "_", deCODE_sexstrata[ss],
                          "_deCode_European_GTh_05062023.txt")
    # Not all pheno-sexstrata combinations have input, so check that the file exists
    if (file.exists(infile_name)) {
      logfile_name <- paste0(main_filepath, "/data/deCode/filtered/", pheno, "_", ss, 
                             ".log")
      outfile_name <- paste0(main_filepath, "/data/deCode/filtered/", pheno, "_", ss, 
                             "_filtered.txt")
      outplot_dir <- paste0(main_filepath, "/plots/deCode/", pheno, "_", ss, "/")
      
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




sexstrata <- c("F", "M", "sex_comb")
phenos <- c("b0", "b1", "k1", "k1_k2", "k1_k2_k3")
strata_plot <- paste0(rep(sexstrata, each = 4), rep(phenos, times = 3))
submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/submit_filter_plot_gwas.sh"

for (ss in strata_plot) {
    job_options <- paste0(
      "--export=",
      paste0(
        "STRATA=\"", ss, "\""
      )
    )
    job_submission <- paste("sbatch", job_options, submission_script)
    system(job_submission)
    print(job_submission)
}

