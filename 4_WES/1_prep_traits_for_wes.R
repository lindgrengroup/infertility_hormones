# Author: Samvida S. Venkatesh
# Date: 21/02/2023

library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

# Read data ----

main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS"

# Covariates file 
covars <- read.table(paste0(main_filepath, "/sample_qc/qcd_sample_file_for_regenie.txt"),
                     sep = "\t", header = T, stringsAsFactors = F)
covars <- covars %>%
  mutate(FID = as.character(FID),
         IID = as.character(IID))

# Hormone files
STRATA <- c("FSH_F", "FSH_M", "FSH_sex_comb",
            "LH_F", "LH_M", "LH_sex_comb", 
            "Oestradiol_F", "Oestradiol_M", "Oestradiol_sex_comb", 
            "Progesterone_F", "Progesterone_sex_comb", 
            "Testosterone_F", "Testosterone_M", "Testosterone_sex_comb")
hr_dat <- lapply(STRATA, function (st) {
  df <- read.table(paste0(main_filepath, "/traits_for_gwas/", st, "_cross_sectional.txt"),
                   sep = "\t", header = T, stringsAsFactors = F)
  colnames(df) <- c("FID", "IID", st)
  df <- df %>%
    mutate(FID = as.character(FID),
           IID = as.character(IID))
  return (df)
})

# Bind columns
hr_dat <- hr_dat %>% reduce(full_join, by = c("FID", "IID"))
full_dat <- full_join(hr_dat, covars, by = c("FID", "IID"))

write.table(full_dat, 
            paste0(main_filepath, "/traits_for_gwas/all_hormones_for_wes_analyses_non_finnish_EUR.txt"),
            sep = "\t", row.names = F, quote = F)

# Plot hormone distributions ----

STRATA <- c("FSH_F", "FSH_M", "FSH_sex_comb",
            "LH_F", "LH_M", "LH_sex_comb", 
            "Oestradiol_F", "Oestradiol_M", "Oestradiol_sex_comb", 
            "Progesterone_F", "Progesterone_sex_comb", 
            "Testosterone_F", "Testosterone_M", "Testosterone_sex_comb")

sex_colour_scheme <- c("#F8766D", "#00BFC4", "#C77CFF")
names(sex_colour_scheme) <- c("F", "M", "sex_comb")

main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS"

full_dat <- read.table(paste0(main_filepath, "/traits_for_gwas/all_hormones_for_wes_analyses_non_finnish_EUR.txt"),
                       sep = "\t", header = T, stringsAsFactors = F)
full_dat$IID <- as.character(full_dat$IID)
id_covars <- read.table(paste0(main_filepath, "/sample_qc/qcd_sample_file_for_regenie.txt"),
                        sep = "\t", header = T, stringsAsFactors = F)
id_covars$IID <- id_covars$IID

getDistPlot <- function (st) {
  plot_dat <- full_dat %>% select(all_of(c("IID", st)))
  plot_dat$sex <- id_covars$sex[match(full_dat$IID,
                                        id_covars$IID)]
  plot_dat <- plot_dat[complete.cases(plot_dat), ]
  
  nf <- length(which(plot_dat$sex == "F"))
  nm <- length(which(plot_dat$sex == "M"))
  
  # Plot the full distribution, regardless of sex
  plot_dat$dummy_sc <- "sex_comb"
  p1 <- ggplot(plot_dat, aes(x = !!as.symbol(st),
                             colour = dummy_sc,
                             fill = dummy_sc)) +
    geom_density(alpha = 0.5) +
    scale_color_manual(values = sex_colour_scheme) +
    scale_fill_manual(values = sex_colour_scheme) +
    labs(x = paste0("RINTed ", st), title = paste0("N = ", nrow(plot_dat)))
  
  p2 <- ggplot(plot_dat, aes(x = !!as.symbol(st),
                             colour = sex,
                             fill = sex)) +
    geom_density(alpha = 0.5) +
    scale_color_manual(values = sex_colour_scheme) +
    scale_fill_manual(values = sex_colour_scheme) +
    labs(x = paste0("RINTed ", st), title = paste0("N_fem = ", nf,
                                                   ", N_male = ", nm))
  
  plot_res <- ggarrange(plotlist = list(p1, p2),
                        ncol = 2, common.legend = T)
  return (plot_res)
}

lapply(STRATA, function (st) {
  p <- getDistPlot(st)
  ggsave(paste0(main_filepath, "/plots/rinted_distribution_", st, ".png"),
         p, units = "in", height = 7, width = 14)
})
