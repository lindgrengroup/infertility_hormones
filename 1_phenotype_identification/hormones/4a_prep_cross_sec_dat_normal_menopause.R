# Author: Samvida S. Venkatesh
# Date: 16/11/2021

library(tidyverse)
theme_set(theme_bw())

# Read data ----

HORMONES <- c("FSH", "LH")
dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/gp_main_data.rds")[HORMONES]
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/230130_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)

# Add in data provider as covariate later if there are enough levels
COVARS <- c("age_event", "age_sq")

# QC log file
qc_log <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/qc/normal_menopause_cross_sectional_qc.txt"

# Remove zero, negative values and get cross-sectional value ----

cross_sec_dat <- lapply(HORMONES, function (hr) {
  # Recover cross-sectional value
  res <- dat[[hr]] %>% group_by(eid) %>%
    filter(value > 0) %>%
    # Get median value and distance of each observation to median
    mutate(medn_value = median(value),
           dist_to_medn = abs(value - median(value))) %>%
    # Arrange by age to get first observation that is closest to median
    group_by(eid) %>% arrange(age_event, .by_group = T) %>%
    slice(which.min(dist_to_medn))
  return (res)
})
names(cross_sec_dat) <- HORMONES

# Add covariates, adjust and RINT values ----

cross_sec_dat <- lapply(HORMONES, function (hr) {
  
  df <- cross_sec_dat[[hr]]
  df$sex <- general_covars$sex[match(df$eid, general_covars$eid)]
  # Only retain individuals that have sex information 
  df <- df[df$sex == "F", ]
  sink(qc_log, append = T)
  cat(paste0("** HORMONE ** ", hr, "\n",
             "\t", "# females: ", nrow(df), "\n"))
  sink()
  
  # Filter to known age at menopause >45yrs
  df$age_at_menopause <- general_covars$age_at_menopause[match(df$eid, 
                                                               general_covars$eid)]
  df <- df[complete.cases(df), ]
  df <- df[df$age_at_menopause > 45, ]
  
  sink(qc_log, append = T)
  cat(paste0("** HORMONE ** ", hr, "\n",
             "\t", "# females with age at menopause >45 years: ", nrow(df), "\n"))
  sink()
  
  # Calculate adjustments and transform (RINT) the value
  # Add covariate for squared age
  df <- df %>% mutate(age_sq = age_event^2)
  # Add data provider as a covariate if there is > 1 provider
  if (length(unique(df$data_provider)) > 1) {
    covars_adj <- c(COVARS, "data_provider")
  }
  
  # Formula for adjustment
  mod_formula <- 
    formula(paste0("value ~ ", paste(covars_adj, collapse = " + ")))
  # Get residuals
  mod_resid <- lm(mod_formula, data = df)$residuals
  
  # RINT
  df$rinted_value <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
  # Write formatted data to file for GWAS if there are enough, i.e. > 1000
  to_write <- data.frame(FID = df$eid, 
                         IID = df$eid,
                         adj_trait = df$rinted_value)
  colnames(to_write) <- c("FID", "IID", hr)
  write.table(to_write,
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/",
                     hr, "_normal_menopause_cross_sectional.txt"),
              sep = "\t", row.names = F, quote = F)
  
})
names(cross_sec_dat) <- HORMONES

saveRDS(cross_sec_dat, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/data/normal_menopause_cross_sec.rds")
