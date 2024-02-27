# Author: Samvida S. Venkatesh
# Date: 17/04/2022

library(tidyverse)

mainpath <- "/well/lindgren/samvida/hormones_infertility/conditional_analysis"

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone") 
SEX_STRATA <- c("F", "M", "sex_comb")
ANCESTRY_GROUPS <- c("EUR", "all_anc")

# Sample sizes ----

sample_sizes <- read.table("/well/lindgren/samvida/hormones_infertility/meta_results_230613/sample_sizes.txt",
                           sep = "\t", header = T, stringsAsFactors = F)

# Lead SNPs ----

lead_snps <- lapply(HORMONES, function (hr) {
  res_list_sex <- lapply(SEX_STRATA, function (sx) {
    res_list_anc <- lapply(ANCESTRY_GROUPS, function (anc_gp) {
      fname <- paste0(mainpath, "/lead_snps/1kg_lead_snp_sumstats_", hr, "_", sx, "_", anc_gp, 
                      "_with_rsids.txt")
      res <- NULL
      if (file.exists(fname)) {
        if (length(readLines(fname)) > 1) {
          res <- read.table(fname, sep = "\t", header = T,
                            stringsAsFactors = F)
          res$ID <- paste0(res$CHROM, "_", res$GENPOS)
          res <- res[, c("ID", "VARID_1KG", "CHROM", "GENPOS")]
          colnames(res) <- c("ID", "VARID_1KG", "chromosome", "position")
        }
      }
      return (res)
    })
    names(res_list_anc) <- ANCESTRY_GROUPS
    return (res_list_anc)
  })
  names(res_list_sex) <- SEX_STRATA
  return (res_list_sex)
})
names(lead_snps) <- HORMONES

# Reported variant lists ----

reported_variants <- lapply(HORMONES, function (hr) {
  res <- read.table(paste0(mainpath, "/published_snps/", hr, "_reported.txt"),
                    sep = " ", header = F, stringsAsFactors = F)
  colnames(res) <- c("ID", "VARID_1KG", "reported_trait")
  res <- res %>%
    mutate(chromosome = gsub(":.*", "", ID),
           position = gsub(".*:", "", ID))
  res$chromosome <- gsub("chr", "", res$chromosome)
  res$chromosome[res$chromosome == "X"] <- 23
  res$ID <- paste0(res$chromosome, "_", res$position)
  return (res[, c("ID", "VARID_1KG", "chromosome", "position")])
})
names(reported_variants) <- HORMONES

reported_variants_all_hormones <- read.table(paste0(mainpath, "/published_snps/all_hormones_reported.txt"),
                                             sep = " ", header = F, stringsAsFactors = F)
colnames(reported_variants_all_hormones) <- c("ID", "VARID_1KG", "reported_trait")
reported_variants_all_hormones <- reported_variants_all_hormones %>%
  mutate(chromosome = gsub(":.*", "", ID),
         position = gsub(".*:", "", ID))
reported_variants_all_hormones$chromosome <- gsub("chr", "", reported_variants_all_hormones$chromosome)
reported_variants_all_hormones$chromosome[reported_variants_all_hormones$chromosome == "X"] <- 23
reported_variants_all_hormones$ID <- paste0(reported_variants_all_hormones$chromosome, "_", reported_variants_all_hormones$position)
reported_variants_all_hormones <- 
  reported_variants_all_hormones[, c("ID", "VARID_1KG", "chromosome", "position")]

# Functions for sorting ----

# Flag and remove any reported variants
flagReportedVars <- function (dat, hormone, log_file) {
  # Unreported for this specific hormone
  unrep_specific <- dat %>% filter(!ID %in% reported_variants[[hormone]]$ID)
  sink(log_file, append = T)
  cat(paste0("\t", "# Previously reported to associate with ", hormone, ": ",
             nrow(dat) - nrow(unrep_specific), 
             "\n"))
  sink()
  
  # Unreported for any hormone
  unrep_any <- dat %>% filter(!ID %in% reported_variants_all_hormones$ID)
  sink(log_file, append = T)
  cat(paste0("\t", "# Previously reported to associate with some hormone: ",
             nrow(dat) - nrow(unrep_any), 
             "\n"))
  sink()
  
  return (list(specific = unrep_specific,
               any = unrep_any))
}

# Calculate +/- 500kb window to perform conditional analysis
WINDOW_SIZE <- 500000 # 500kb window
getCondList <- function (hormone, variant_chr, variant_pos) {
  if (hormone == "any") {
    condSNPS <- reported_variants_all_hormones %>% filter(chromosome == variant_chr &
                                                            position >= variant_pos - WINDOW_SIZE &
                                                            position <= variant_pos + WINDOW_SIZE)
  } else {
    condSNPS <- reported_variants[[hormone]] %>% filter(chromosome == variant_chr &
                                                          position >= variant_pos - WINDOW_SIZE &
                                                          position <= variant_pos + WINDOW_SIZE)
  }
  if (nrow(condSNPS) > 0) {
    condSNPS <- condSNPS$VARID_1KG
  } else condSNPS <- NA
  return (condSNPS)
}

# Write novel SNPs and published SNPs in window
writeNovelPublished <- function (unrep_df, hormone, path_store) {
  # Specifically unreported for this hormone
  novel_snps <- c()
  for (i in 1:nrow(unrep_df)) {
    pub_snps_in_window <- getCondList(hormone = hormone,
                                      unrep_df$chromosome[i],
                                      unrep_df$position[i])
    # If there is nothing in the window, report this as a novel SNP
    if (all(is.na(pub_snps_in_window))) {
      novel_snps <- c(novel_snps, unrep_df$ID[i])
    } else {
      # Add SNP itself to the list
      pub_snps_in_window <- unique(c(unrep_df$VARID_1KG[i], 
                                     pub_snps_in_window))
      write.table(pub_snps_in_window,
                  paste0(path_store, "/", hormone, "_published_snps_around_",
                         unrep_df$VARID_1KG[i], ".txt"),
                  sep = "\t", row.names = F, col.names = F, quote = F)
    }
  }
  return (novel_snps)
}

SUBMISSION_SCRIPT <- "/well/lindgren/samvida/hormones_infertility/scripts/gcta_cojo_conditional_analysis.sh"
# Submit GCTA-COJO conditional analysis jobs
submitGCTA <- function (hr, sx, anc, 
                        unrep_df, hormone_type) {
  snps_write <- unrep_df$VARID_1KG
  write.table(snps_write,
              paste0(mainpath, "/unreported_snps/", hr, "_", sx, "_",
                     anc, "/unreported_", hormone_type, ".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  job_options <- paste0(
    "--export=",
    paste0(
      "HORMONE=\"", hr, "\",",
      "SEX_STRATA=\"", sx, "\",",
      "ANC_GROUP=\"", anc, "\",",
      "SAMPLE_SIZE=\"", sample_sizes$N[sample_sizes$strata == paste0(hr, "_", sx, "_", anc)], "\",",
      "HORMONE_TYPE=\"", hormone_type, "\""
    )
  )
  job_submission <- paste("sbatch", job_options, SUBMISSION_SCRIPT)
  system(job_submission)
  print(job_submission)
}

# Apply across all strata ----

lapply(HORMONES, function (hr) {
  lapply(SEX_STRATA, function (sx) {
    lapply(ANCESTRY_GROUPS, function (anc_gp) {
      df <- lead_snps[[hr]][[sx]][[anc_gp]]
      # If the strata exists
      if (!is.null(df)) {
        log_file <- paste0(mainpath, "/logs/", hr, "_", sx, "_", anc_gp, ".txt")
        tmp_dir <- paste0(mainpath, "/unreported_snps/", hr, "_", sx, "_", anc_gp)
        dir.create(tmp_dir)
        
        # Record number of SNPs
        sink(log_file, append = T)
        cat(paste0("# Lead SNPs: ", nrow(df), "\n"))
        sink()
        
        # Flag and remove reported SNPs
        unreported_df <- flagReportedVars(dat = df, hormone = hr, log_file = log_file)
        
        # If there are any SNPs left, get SNPs to perform conditional analysis
        
        # For this specific hormone
        if (nrow(unreported_df[["specific"]]) > 0) {
          novel_specific <- writeNovelPublished(unreported_df[["specific"]],
                                                hormone = hr,
                                                path_store = tmp_dir)
          write.table(novel_specific,
                      paste0(mainpath, "/novel_snps/", hr, "_", sx, "_", anc_gp, 
                             "_novel_specific.txt"),
                      sep = "\t", row.names = F, col.names = F, quote = F)
          
          # Submit GCTA-COJO jobs
          # Remove novel SNPs from the dataframes first
          poss_novel_specific <- unreported_df[["specific"]]
          poss_novel_specific <- poss_novel_specific %>% filter(!ID %in% novel_specific)
          
          if (nrow(poss_novel_specific) > 0) {
            submitGCTA(hr, sx, anc_gp,
                       poss_novel_specific, hormone_type = hr)
          }
        }
        
        # For any hormone
        if (nrow(unreported_df[["any"]]) > 0) {
          novel_any <- writeNovelPublished(unreported_df[["any"]],
                                           hormone = "any",
                                           path_store = tmp_dir)
          write.table(novel_any,
                      paste0(mainpath, "/novel_snps/", hr, "_", sx, "_", anc_gp, 
                             "_novel_any.txt"),
                      sep = "\t", row.names = F, col.names = F, quote = F)
          
          # Submit GCTA-COJO jobs
          # Remove novel SNPs from the dataframes first
          poss_novel_any <- unreported_df[["any"]]
          poss_novel_any <- poss_novel_any %>% filter(!ID %in% novel_any)
          
          if (nrow(poss_novel_any) > 0) {
            submitGCTA(hr, sx, anc_gp,
                       poss_novel_any, hormone_type = "any")
          }
        }
      }
    })
  })
})
