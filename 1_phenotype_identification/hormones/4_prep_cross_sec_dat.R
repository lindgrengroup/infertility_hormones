# Author: Samvida S. Venkatesh
# Date: 16/11/2021

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

# Read data ----

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone")
dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/gp_main_data.rds")[HORMONES]
non_wb_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/non_wb_gp_main_data.rds")[HORMONES]
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/230130_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)

nfe_ids <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_non_finnish_eur.txt",
                      sep = "\t", header = F, stringsAsFactors = F)$V1
nfe_ids <- nfe_ids[nfe_ids > 0]
nfe_ids <- as.character(nfe_ids)

SEX_STRATA <- c("F", "M", "sex_comb")
# Add in data provider as covariate later if there are enough levels
COVARS <- c("age_event", "age_sq")

# QC log file
qc_log <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/qc/full_nfe_cross_sectional_qc.txt"

# Remove zero, negative values and get cross-sectional value ----

cross_sec_dat <- lapply(HORMONES, function (hr) {
  full_dat <- bind_rows(dat[[hr]], non_wb_dat[[hr]]) %>%
    select(-year_of_birth) %>%
    filter(eid %in% nfe_ids)
  
  # Recover cross-sectional value
  res <- full_dat %>% group_by(eid) %>%
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

# Plot the distribution of each hormone value ----

colpalette_use <- c("#F8766D", "#619CFF")
names(colpalette_use) <- c("F", "M")

meno_colpalette <- c("#40004B", "#00441B")
names(meno_colpalette) <- c("post_menopause", "pre_menopause")

lapply(HORMONES, function (hr) {
  for_plot <- cross_sec_dat[[hr]] 
  for_plot$sex <- general_covars$sex[match(for_plot$eid,
                                           general_covars$eid)]
  lapply(SEX_STRATA, function (sx) {
    if (sx != "sex_comb") for_plot <- for_plot %>% filter(sex == sx)
    
    for_plot <- for_plot %>% filter(!is.na(sex))
    
    p1 <- ggplot(for_plot,
                 aes(x = value, fill = sex, colour = sex)) +
      geom_density(alpha = 0.2) +
      scale_color_manual(values = colpalette_use) +
      scale_fill_manual(values = colpalette_use) + 
      labs(x = paste0(hr))
    
    p2 <- ggplot(for_plot,
                 aes(x = value, fill = sex, colour = sex)) +
      geom_density(alpha = 0.2) +
      scale_color_manual(values = colpalette_use) +
      scale_fill_manual(values = colpalette_use) + 
      scale_x_continuous(trans = "log") +
      labs(x = paste0(hr, ", log axis"))
    
    resplot <- ggarrange(p1, p2, ncol = 2)
    ggsave(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/plots/", 
                  hr, "_", sx, ".png"), 
           resplot)
    
  })
})

# Add covariates, adjust and RINT values ----

cross_sec_dat <- lapply(HORMONES, function (hr) {
  
  df <- cross_sec_dat[[hr]]
  df$sex <- general_covars$sex[match(df$eid, general_covars$eid)]
  # Only retain individuals that have sex information 
  df <- df[complete.cases(df), ]
  sink(qc_log, append = T)
  cat(paste0("** HORMONE ** ", hr, "\n",
             "\t", "# males: ", sum(df$sex == "M"), "\n",
             "\t", "# females: ", sum(df$sex == "F"), "\n"))
  sink()
  
  # Calculate adjustments and transform (RINT) the value
  res <- lapply(SEX_STRATA, function (sx) {
    if (sx == "sex_comb") {
      sub_df <- df
      covars_adj <- c(COVARS, "sex")
    } else {
      sub_df <- df %>% filter(sex == sx)
      covars_adj <- COVARS
    }
    
    # Add covariate for squared age
    sub_df <- sub_df %>% mutate(age_sq = age_event^2)
    # Add data provider as a covariate if there is > 1 provider
    if (length(unique(sub_df$data_provider)) > 1) {
      covars_adj <- c(covars_adj, "data_provider")
    }
    
    # Formula for adjustment
    mod_formula <- 
      formula(paste0("value ~ ", paste(covars_adj, collapse = " + ")))
    # Get residuals
    mod_resid <- lm(mod_formula, data = sub_df)$residuals
    
    # Plot by menopause status again to see if the adjustment 
    # has taken care of the peak differences
    
    if (sx == "F") {
      sub_df$adj_value <- mod_resid
      sub_df$age_at_menopause <- general_covars$age_at_menopause[match(sub_df$eid,
                                                                       general_covars$eid)]
      
      menopause_plot <- sub_df %>% 
        filter(!is.na(age_at_menopause)) %>%
        mutate(menopause_status = ifelse(age_event <= age_at_menopause,
                                         "pre_menopause", "post_menopause"))
      
      p1_meno <- ggplot(menopause_plot,
                          aes(x = value, fill = menopause_status, 
                              colour = menopause_status)) +
        geom_density(alpha = 0.2) +
        scale_color_manual(values = meno_colpalette) +
        scale_fill_manual(values = meno_colpalette) + 
        scale_x_continuous(trans = "log") +
        labs(x = paste0("Unadjusted ", hr, ", log axis"))
      
      p2_meno <- ggplot(menopause_plot,
                        aes(x = adj_value, fill = menopause_status, 
                            colour = menopause_status)) +
        geom_density(alpha = 0.2) +
        scale_color_manual(values = meno_colpalette) +
        scale_fill_manual(values = meno_colpalette) + 
        scale_x_continuous(trans = "log") +
        labs(x = paste0("Age and age-sq adjusted ", hr, ", log axis"))
      
      meno_plot <- ggarrange(p1_meno, p2_meno, ncol = 2)
      
      ggsave(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/plots/", 
                    hr, "_by_menopause_status.png"), 
             meno_plot)
      
    }
    
    # RINT
    sub_df$rinted_value <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
    # Write formatted data to file for GWAS if there are enough, i.e. > 1000
    if (nrow(sub_df) > 1000) {
      to_write <- data.frame(FID = sub_df$eid, 
                             IID = sub_df$eid,
                             adj_trait = sub_df$rinted_value)
      colnames(to_write) <- c("FID", "IID", hr)
      write.table(to_write,
                  paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/",
                         hr, "_", sx, "_cross_sectional.txt"),
                  sep = "\t", row.names = F, quote = F)
    } else {
      to_write <- NULL
    }
    
    return (to_write)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(cross_sec_dat) <- HORMONES

saveRDS(cross_sec_dat, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/data/split_sex_cross_sec.rds")
