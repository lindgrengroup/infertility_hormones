# Author: Samvida S. Venkatesh
# Date: 31/08/21

library(tidyverse)
library(lubridate)

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
ukb_path <- "" # REDACTED

# Log file ----

indiv_qc_log <- paste0(mainpath, "/qc/gp_and_main_longit_filter.txt")

# Read data ----

# WB ancestry ids
wb_ids <- read.table(paste0(gen_resources_path, "/eids_white_british.txt"),
                     sep = "\t", header = F, stringsAsFactors = F)$V1

# UK Biobank main phenotype file
pheno <- read.table(paste0(ukb_path, "/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv"),
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "eid"
# Retain only WB ancestry
pheno <- subset(pheno, pheno$eid %in% wb_ids)

# Biomarkers file from UKBB
biomarkers <- read.table(paste0(ukb_path, "/DATA/Biomarker_data/ukb27722.csv"),
                         header = T, sep = ",", na.string = c("NA", "", "."), 
                         stringsAsFactors = F)
colnames(biomarkers) <- gsub("X", "f.", colnames(biomarkers))
colnames(biomarkers)[1] <- "eid"
# Retain only WB ancestry
biomarkers <- subset(biomarkers, biomarkers$eid %in% wb_ids)

# Biomarker data from primary care
gp_dat <- readRDS(paste0(mainpath, "/data_passed_popn_QC.rds"))
PHENOTYPES <- names(gp_dat)

# Main phenotype or biomarker file field codes
field_codes <- read.table(paste0(mainpath, "/code_lists/main_phenotypes/ukb_main_biomarker_codelist.txt"),
                          sep = "\t", header = T, stringsAsFactors = F)
MAIN_PHENOS <- field_codes$biomarker[field_codes$file == "Main"]
BIOM_PHENOS <- field_codes$biomarker[field_codes$file == "Biomarkers"]
ONLY_GP <- PHENOTYPES[!PHENOTYPES %in% c(MAIN_PHENOS, BIOM_PHENOS)]

# Add relevant UKB data to each phenotype ----

# Function to convert main pheno wide files to long format
wide_to_long <- function (wide_df, p) {
  res <- wide_df
  # Merge in centre codes 
  extra_codes <- pheno[, c("eid", 
                           "f.54.0.0", "f.54.1.0", "f.54.2.0")]
  res <- merge(res, extra_codes, by = "eid")
  # Pivot to long format
  # Special case for WHR
  if (p == "WHR") {
    colnames(res) <- c("eid", 
                       "event_dt.1", "event_dt.2", "event_dt.3",
                       "waist.1", "waist.2", "waist.3", 
                       "hip.1", "hip.2", "hip.3",
                       "data_provider.1", "data_provider.2", 
                       "data_provider.3")
    res <- pivot_longer(res, 
                        cols = -1, 
                        names_to = c(".value", "set"),
                        names_pattern = "(.+).(.+)", 
                        values_drop_na = T) %>%
      mutate(value = waist/hip)
  } else if (p %in% c("FEV1", "FVC")) {
    # Match column names with GP data, FEV1 and FVC only have 1 value
    colnames(res) <- c("eid", "event_dt", "value", 
                       "data_provider")
  } else if (p %in% MAIN_PHENOS)  {
    # Match column names with GP data, main phenos have 2 repeats
    colnames(res) <- c("eid", 
                       "event_dt.1", "event_dt.2", "event_dt.3",
                       "value.1", "value.2", "value.3",
                       "data_provider.1", "data_provider.2", 
                       "data_provider.3")
    res <- pivot_longer(res, cols = -1, 
                        names_to = c(".value", "set"),
                        names_pattern = "(.+).(.+)", 
                        values_drop_na = T)
    
  } else if (p %in% BIOM_PHENOS) {
    # Match column names with GP data, biom phenos have 1 repeat
    colnames(res) <- c("eid", 
                       "event_dt.1", "event_dt.2", 
                       "value.1", "value.2", 
                       "data_provider.1", "data_provider.2", 
                       "data_provider.3")
    res <- pivot_longer(res, cols = -1, 
                        names_to = c(".value", "set"),
                        names_pattern = "(.+).(.+)", 
                        values_drop_na = T)
  }
  
  res$data_provider <- paste0("UKBB", res$data_provider) 
  res$biomarker <- p
  res <- res[, c("eid", "data_provider", "event_dt",
                 "value", "biomarker")]
  return (res)
}

# Function to get age at event from birth date
get_age <- function (long_df) {
  # Get date of birth from main phenotype file (month and year)
  tmp <- merge(long_df, pheno[, c("eid", "f.34.0.0", "f.52.0.0")],
               by = "eid")
  colnames(tmp) <- c(colnames(long_df), "year_of_birth", "month_of_birth")
  # Assign date of birth (first of the month) based on month and year of birth 
  tmp$dob <- as.Date(paste(tmp$year_of_birth, tmp$month_of_birth,
                           1, sep = "-"))
  
  # Calculate age at event
  tmp$event_dt <- as.Date(tmp$event_dt, "%Y-%m-%d")
  tmp$dob <- as.Date(tmp$dob, "%Y-%m-%d")
  # Remove inconsistencies (NA in event date, event date before date-of-birth, 
  # and event date after 2020)
  tmp <- tmp %>% filter(!is.na(event_dt) & event_dt > dob & 
                           event_dt <= as.Date("2020-10-01"))
  tmp$age_event <- interval(tmp$dob, tmp$event_dt) / years(1)
  
  res <- tmp[, c("eid", "data_provider", "event_dt", "age_event",
                 "value", "biomarker")]
  return (res)
}

# Function to standardise column types 
standardise_coltypes <- function (long_df) {
  res <- long_df[, c("eid", "data_provider", 
                     "event_dt", "age_event", "value",
                     "biomarker")]
  # Convert all columns to standard type
  res$eid <- factor(long_df$eid)
  res$data_provider <- factor(long_df$data_provider)
  res$event_dt <- as.Date(long_df$event_dt, "%Y-%m-%d")
  res$age_event <- as.numeric(long_df$age_event)
  res$value <- as.numeric(long_df$value)
  return (res)
}

gp_and_main <- lapply(PHENOTYPES, function (p) {
  gp_df <- gp_dat[[p]]
  # Only keep individuals with at least one trait measurement in GP
  ids_keep <- unique(gp_df$eid)
  print(paste0("Running phenotype: ", p))
  
  # If the phenotype isn't in any UKB files
  if (p %in% ONLY_GP) {
    res <- gp_df
  } else {
    # Check whether phenotype is in main, biomarkers, or WHR
    # Special case for WHR as it has to be calculated from WC and HC
    if (p == "WHR") {
      value_codes <- c("f.48.0.0", "f.48.1.0", "f.48.2.0",
                       "f.49.0.0", "f.49.1.0", "f.49.2.0")
      date_codes <- c("f.53.0.0", "f.53.1.0", "f.53.2.0")
      main_df <- pheno[, c("eid", date_codes, value_codes)]
    } else if (p %in% c("FEV1", "FVC")) {
      # Get relevant columns for value, date, and assessment centre
      value_codes <- 
        colnames(biomarkers)[grep(paste0("f.", 
                                         field_codes$field[field_codes$biomarker == p]),
                                  colnames(biomarkers))]
      # Get date column from main phenotype file, only first date because 
      # there is only 1 value for FEV1 and FVC
      date_codes <- c("f.53.0.0")
      main_df <- pheno[, c("eid", date_codes)]
      main_df <- merge(main_df, biomarkers[, c("eid", value_codes)],
                       by = "eid")
    } else if (p %in% MAIN_PHENOS) {
      # Get relevant columns for value, date, and assessment centre
      value_codes <- 
        colnames(pheno)[grep(paste0("f.", 
                                    field_codes$field[field_codes$biomarker == p]),
                             colnames(pheno))]
      date_codes <- c("f.53.0.0", "f.53.1.0", "f.53.2.0")
      main_df <- pheno[, c("eid", date_codes, value_codes)]
    } else if (p %in% BIOM_PHENOS) {
      value_codes <- 
        colnames(biomarkers)[grep(paste0("f.", 
                                         field_codes$field[field_codes$biomarker == p]),
                                  colnames(biomarkers))]
      date_codes <- 
        colnames(biomarkers)[grep(paste0("f.", field_codes$date_field[field_codes$biomarker == p]),
                                  colnames(biomarkers))]
      main_df <- biomarkers[, c("eid", date_codes, value_codes)]
    }
    # Convert main df to long format 
    main_df <- subset(main_df, main_df$eid %in% ids_keep)
    main_df <- wide_to_long(main_df, p)
    # Get age at event from birth date
    main_df <- get_age(main_df)
    # Standardise columns
    main_df <- standardise_coltypes(main_df)
    # Combine with GP data
    res <- bind_rows(gp_df, main_df) %>% distinct()
  }
  return (res)
})
names(gp_and_main) <- PHENOTYPES

saveRDS(gp_and_main, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/gp_main_data_passed_longit_filter.rds")

