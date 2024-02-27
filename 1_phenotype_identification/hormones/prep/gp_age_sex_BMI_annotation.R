# Author: Samvida S. Venkatesh
# Date: 04/12/20

# Create new UKBB GP clinical file with age at each event and mean UKBIOBANK BMI

# NEEDS AT LEAST 5 CORES (PREFERABLY 7 CORES)

NPAR_CORES <- 5

library(lubridate)
library(eeptools)
library(foreach)
library(doParallel) 
library(dplyr)

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED

# Read data ----

# Read gp_clinical file
gp_clinical <- read.table(paste0(mainpath, "/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt"),
                          sep = "\t", header = T, comment.char = "$",
                          stringsAsFactors = F)

# Remove individuals who have withdrawn consent
withdrawn <- read.table(paste0(mainpath, "/QC/w11867_20200820.csv"), 
                        header = F)
gp_clinical <- subset(gp_clinical, !gp_clinical$eid %in% withdrawn$V1)

# Main UKBB phenotypes file
pheno <- read.table(paste0(mainpath, "/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv"),
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# Add info from main UKBB phenotypes file ----

## Add information on date-of-birth and sex ----

gp_clinical <- merge(gp_clinical, pheno[, c("f.eid", "f.31.0.0", 
                                            "f.34.0.0", "f.52.0.0")], 
                     by.x = "eid", by.y = "f.eid", all.x = T)
colnames(gp_clinical)[9:11] <- c("sex", "year_of_birth", "month_of_birth")

# Assign date of birth (first of the month) based on month and year of birth 
gp_clinical$dob <- as.Date(paste(gp_clinical$year_of_birth, 
                                 gp_clinical$month_of_birth,
                                 1, sep = "-"))
# Change sex coding to M, F
gp_clinical$sex <- ifelse(gp_clinical$sex == 0, "F", "M")

## Add mean UKBIOBANK BMI ----

ukbb_bmi <- pheno[, c("f.eid", 
                      "f.21001.0.0", "f.21001.1.0", "f.21001.2.0")]
ukbb_bmi$mean_UKBB_BMI <- rowMeans(ukbb_bmi[, 2:4], na.rm = T)
colnames(ukbb_bmi)[1] <- "eid"

gp_clinical <- merge(gp_clinical, ukbb_bmi[, c("eid", "mean_UKBB_BMI")],
                 by = "eid", all.x = T)

# Calculate individual age at event ----

gp_clinical$event_dt <- as.Date(gp_clinical$event_dt, "%Y-%m-%d")
gp_clinical$dob <- as.Date(gp_clinical$dob, "%Y-%m-%d")

# Remove inconsistencies (NA in event date, event date before date-of-birth, 
# and event date after 2020)
cleaned <- subset(gp_clinical, !is.na(gp_clinical$event_dt))
cleaned <- subset(cleaned, cleaned$event_dt > cleaned$dob)
cleaned <- subset(cleaned, cleaned$event_dt <= as.Date("2020-10-01"))

# Chunk this into smaller pieces for computational
# efficiency; store smaller chunks and delete later

CSIZE <- 1000000 # chunk size 1 million
n <- nrow(cleaned)
r <- rep(1:ceiling(n / CSIZE), each = CSIZE, length.out = n)[1:n]
d <- split(cleaned, r)

# Parallel processing

cl <- makeCluster(NPAR_CORES)
registerDoParallel(cl)

foreach(i = 1:length(d), .packages = c("eeptools")) %dopar% {
  
  # Age calculation
  p <- d[[i]]
  p$age_years <- age_calc(p$dob, p$event_dt, units = "years")
  
  # Keep track of count
  iter <- paste("agecalc_p", i, sep = "")         
  
  # Save (temporary)
  saveRDS(p, paste0(gen_resources_path, "/temp_files/", 
                   iter, ".rds"))
  
  # Forget whatever was in memory since we've saved the RDS
  return(1) 
}

stopCluster(cl)

# Rebuild full GP clinical file

iterations <- list.files(paste0(gen_resources_path, "/temp_files"))

pdata <- lapply(iterations, function(i) {
  readRDS(paste0(gen_resources_path, "/temp_files/", i))
})

pdata <- bind_rows(pdata)

write.table(pdata, paste0(gen_resources_path, "/gp_clinical_annotated.txt"),
            sep = "\t", quote = F, row.names = F)

