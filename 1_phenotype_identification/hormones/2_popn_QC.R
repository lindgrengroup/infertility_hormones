# Author: Samvida S. Venkatesh
# Date: 30/08/21

library(tidyverse)

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED

# Read and subset required data ----

# All quantitative traits
dat <- readRDS(paste0(mainpath, "/data_passed_primary_care_n_qc.rds"))
PHENOTYPES <- names(dat)

# Bariatric surgery codes
bariatric_codes <- read.table(paste0(gen_resources_path, "/bariatric_surgery_records.txt"),
                              sep = "\t", header = T, comment.char = "$",
                              stringsAsFactors = F)
bariatric_codes$event_dt <- as.Date(bariatric_codes$event_dt, "%Y-%m-%d")

# Implausible values from UKB min/max
ukb_ranges <- read.table(paste0(mainpath, "/qc/ukb_min_max.txt"),
                         sep = "\t", header = T, 
                         stringsAsFactors = F)

# Functions to clean data ----

## Age ----

ageFilter <- function (df, qc_log_file) {
  
  remove_young <- df$eid[df$age_event < 20]
  remove_old <- df$eid[df$age_event > 80]
  
  cleaned <- subset(df, df$age_event >= 20)
  cleaned <- subset(cleaned, cleaned$age_event <= 80)
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Age at measurement < 20 years: ", "\n",
             "\t", "Number of measurements = ", 
             length(remove_young), "\n",
             "\t", "Number of individuals = ", 
             length(unique(remove_young)), "\n",
             "**FILTER** EXCLUDED, Age at measurement > 80 years: ", "\n",
             "\t", "Number of measurements = ", 
             length(remove_old), "\n",
             "\t", "Number of individuals = ", 
             length(unique(remove_old)), "\n"))
  sink()
  
  return (cleaned)
}

## Bariatric surgery ----

barSurgeryFilter <- function (df, qc_log_file) {
  
  # Exclude individuals with history of bariatric surgery
  bar_hist_EIDS <- bariatric_codes$eid[is.na(bariatric_codes$event_dt)]
  # Record date of surgery for individuals with that information
  df$post_bariatric_surgery_flag <- 
    bariatric_codes$event_dt[match(df$eid, bariatric_codes$eid)]
  df$post_bariatric_surgery_flag <- df$event_dt > 
    df$post_bariatric_surgery_flag | df$eid %in% bar_hist_EIDS
  df$post_bariatric_surgery_flag[is.na(df$post_bariatric_surgery_flag)] <- F
  
  cleaned <- subset(df, !df$post_bariatric_surgery_flag)
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Measured post-bariatric surgery: ", "\n",
             "\t", "Number of measurements = ", 
             sum(df$post_bariatric_surgery_flag, na.rm = T), "\n",
             "\t", "Number of individuals = ", 
             length(unique(df$eid)) - length(unique(cleaned$eid)), "\n"))
  sink()
  
  return (cleaned)
  
}

## Implausible values ----

implausibleFilter <- function (df, p, qc_log_file) {
  # Remove implausible values defined by +/-10% of UKB
  IMP_MIN <- ukb_ranges$min[ukb_ranges$biomarker == p]
  IMP_MAX <- ukb_ranges$max[ukb_ranges$biomarker == p]
  
  if (length(IMP_MIN) > 0 & length(IMP_MAX) > 0) {
    remove_mins <- df[df$value < IMP_MIN, ]
    remove_maxes <- df[df$value > IMP_MAX, ]
    
    cleaned <- subset(df, df$value >= IMP_MIN & 
                        df$value <= IMP_MAX)
    
    # Report QC metrics
    sink(qc_log_file, append = T)
    cat(paste0("**FILTER** EXCLUDED, Implausible value below ", 
               IMP_MIN, ":", "\n",
               "\t", "Number of measurements = ", 
               dim(remove_mins)[1], "\n",
               "\t", "Number of individuals = ", 
               sum(!unique(remove_mins$eid) %in% unique(cleaned$eid)), "\n",
               "**FILTER** EXCLUDED, Implausible value above ", 
               IMP_MAX, ":", "\n",
               "\t", "Number of measurements = ", 
               dim(remove_maxes)[1], "\n",
               "\t", "Number of individuals = ", 
               sum(!unique(remove_maxes$eid) %in% unique(cleaned$eid)), "\n"))
    sink()
  } else {
    cleaned <- df
    sink(qc_log_file, append = T)
    cat(paste0("**No UKB implausible values**", "\n"))
    sink()
  }
  return (cleaned)
}

## Extreme values ----

extremeFilter <- function (df, qc_log_file) {
  # Remove values +/- 5 S.D. away from the population mean
  popn_mean <- mean(df$value, na.rm = T)
  popn_sd <- sd(df$value, na.rm = T)
  max_outlier <- popn_mean + 5*popn_sd
  min_outlier <- popn_mean - 5*popn_sd
  cleaned <- subset(df, df$value <= max_outlier & df$value >= min_outlier)
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Value > 5 S.D. away from population mean: ", 
             "\n",
             "\t", "Number of measurements = ", 
             dim(df)[1] - dim(cleaned)[1], "\n",
             "\t", "Number of individuals = ", 
             length(unique(df$eid)) - length(unique(cleaned$eid)), "\n"))
  sink()
  return (cleaned)
}


# Apply QC filters ----

cleaned_dat <- lapply(PHENOTYPES, function (p) {
  
  log_file_p <- paste0(mainpath, "/qc/popn_qc/log_file_",
                       p, ".txt")
  cleaned <- ageFilter(dat[[p]], log_file_p)
  if (p %in% c("BMI", "Weight", "WC", "WHR")) {
    cleaned <- barSurgeryFilter(cleaned, log_file_p)
  }
  cleaned <- implausibleFilter(cleaned, p, log_file_p)
  cleaned <- extremeFilter(cleaned, log_file_p)
  
  return (cleaned)
  
})
names(cleaned_dat) <- PHENOTYPES

saveRDS(cleaned_dat, paste0(mainpath, "/data_passed_popn_QC.rds"))
