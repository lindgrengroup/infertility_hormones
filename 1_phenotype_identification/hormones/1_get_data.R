# Author: Samvida S. Venkatesh
# Date: 27/08/21

library(tidyverse)
library(lubridate)

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED

# Log file ----

qc_log <- paste0(mainpath, "/qc/num_repeat_measures.txt")

# Read data ----

# List of biomarkers
BM_LIST <- read.table(paste0(mainpath, "/code_lists/traits_available.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)$biomarker
  
# WB ancestry ids
wb_ids <- read.table(paste0(gen_resources_path, "/eids_white_british.txt"),
                       sep = "\t", header = F, stringsAsFactors = F)$V1

# Annotated GP data
gp_clinical <- read.table(paste0(gen_resources_path, "/gp_clinical_annotated.txt"),
                          sep = "\t", header = T, comment.char = "$",
                          stringsAsFactors = F)

# Retain WB subset
gp_clinical <- subset(gp_clinical, gp_clinical$eid %in% wb_ids)

# Extract relevant read codes ----

bm_gp <- lapply(BM_LIST, function (bm) {
  # Get read codes for biomarker
  codes <- read.table(paste0(mainpath, "/code_lists/primary_care/",
                             bm, "_cleaned_refined_final.csv"), 
                      sep = ",", header = T, stringsAsFactors = F)
  
  R2codes <- codes$readcode[codes$terminology == "read2"]
  R3codes <- codes$readcode[codes$terminology == "ctv3"]
  
  # Create regexp to get codes
  vrep <- Vectorize(rep.int, "times", SIMPLIFY = T)
  
  ndiff <- 5 - nchar(R2codes)
  add_ons <- vrep(".", times = ndiff)
  add_ons <- lapply(add_ons, function (x) paste0(x, collapse = ""))
  R2codes <- paste0(R2codes, unlist(add_ons))
  
  ndiff <- 5 - nchar(R3codes)
  add_ons <- vrep(".", times = ndiff)
  add_ons <- lapply(add_ons, function (x) paste0(x, collapse = ""))
  R3codes <- paste0(R3codes, unlist(add_ons))
  
  # Subset relevant primary care data
  if (length(R2codes) > 0) {
    R2grep <- paste(R2codes, collapse = "|")
    r2_rows <- grep(R2grep, gp_clinical$read_2)
  } else { r2_rows <- numeric(0) }
  
  if (length(R3codes) > 0) {
    R3grep <- paste(R3codes, collapse = "|")
    r3_rows <- grep(R3grep, gp_clinical$read_3)
  } else { r3_rows <- numeric(0) }
  
  all_rows <- unique(c(r2_rows, r3_rows))
  
  if (length(all_rows) == 0) {
    print(paste0("No values for biomarker: ", bm))
    df <- data.frame(eid = character(0), data_provider = character(0), 
                     event_dt = character(0), age_event = numeric(0), 
                     value = numeric(0), biomarker = character(0))
    
  } else {
    # Retain only the columns of interest
    # Only retain value1 as value2 and value3 contain measurement units
    # or misleading values
    df <- gp_clinical[all_rows, 
                      c("eid", "data_provider", "event_dt", "age_years", "value1")]
    colnames(df) <- c("eid", "data_provider", "event_dt", "age_event", "value")
    
    # Remove negative and non-numeric values in value column
    df$value <- as.numeric(df$value)
    df <- subset(df, !is.na(df$value))
    df <- subset(df, df$value >= 0)
    if (nrow(df) == 0) {
      print(paste0("No non-zero values for biomarker: ", bm))
      df <- data.frame(eid = character(0), data_provider = character(0), 
                       event_dt = character(0), age_event = numeric(0), 
                       value = numeric(0), biomarker = character(0))
    } else {
      df$data_provider <- paste0("GP", df$data_provider)
      df$biomarker <- bm
    }
  }
  # Convert all columns to standard type
  df$eid <- factor(df$eid)
  df$data_provider <- factor(df$data_provider)
  df$event_dt <- as.Date(df$event_dt, "%Y-%m-%d")
  df$age_event <- as.numeric(df$age_event)
  df$value <- as.numeric(df$value)
  return (df)
})

names(bm_gp) <- BM_LIST

# Initial QC ----

# Retain only biomarkers with > 1 measurement in > 1000 
# (or 100 for hormones + anthro) individuals in primary care 

HORMONES <- c("17_OHP", "AMH", "DHEAS", "DHT", "FAI", "FSH", "GnRH",
              "HCG", "HPL", "Inhibin", "LH", "Oestradiol", "Oestriol",
              "Oestrone", "Osteocalcin", "Oxytocin", "Progesterone",
              "Prolactin", "Relaxin", "Testosterone",
              "BMI", "WC", "WHR", "Weight")

summarise_nmeasures <- bind_rows(bm_gp)

summarise_repeat <- lapply(BM_LIST, function (bm) {
  # Count number of repeat measures per individual
  df <- bm_gp[[bm]]
  summ_res <- df %>% group_by(eid) %>% summarise(nmeasures = n())
  ct <- length(which(summ_res$nmeasures > 1))
  res <- data.frame(biomarker = bm, 
                    nindivs_with_repeat = ct)
  return (res)
})
summarise_repeat <- bind_rows(summarise_repeat)
# Write table of number of measures 
write.table(summarise_repeat, qc_log,
            sep = "\t", row.names = F, quote = F)

# Biomarkers to retain:
summarise_repeat$retain <- F
summarise_repeat$retain <- 
  ifelse(summarise_repeat$biomarker %in% HORMONES, 
         summarise_repeat$nindivs_with_repeat > 100,
         summarise_repeat$nindivs_with_repeat > 1000)

RETAIN <- summarise_repeat$biomarker[summarise_repeat$retain]

bm_gp <- bm_gp[RETAIN]

saveRDS(bm_gp, 
        paste0(mainpath, "/data_passed_primary_care_n_qc.rds"))

