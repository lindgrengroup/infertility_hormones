# Author: Samvida S. Venkatesh
# Date: 16/11/2021

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

# QC log file
qc_log <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/qc/fsh_lh_age_menopause.txt"

# Read data ----

HORMONES <- c("FSH", "LH")
hormone_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/data/split_sex_cross_sec.rds")[HORMONES]
# Only keep female data
hormone_dat <- lapply(hormone_dat, function (df) {
  df[["F"]]
})

# Get age at menopause
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/230130_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)
aam_dat <- general_covars %>%
  select(all_of(c("eid", "age_at_menopause"))) %>%
  filter(!is.na(age_at_menopause))

sink(qc_log, append = T)
cat(paste0("# female participants with age at menopause: ", 
           nrow(aam_dat), "\n"))
cat(paste0("\t", 
           "# mean (SE) age at menopause: ", 
           mean(aam_dat$age_at_menopause), " (",  
           sd(aam_dat$age_at_menopause)/sqrt(nrow(aam_dat)), ")", "\n"))
sink()

# Wrangle data to add presence of hormone ----

dat_plot <- aam_dat %>%
  mutate(has_FSH = eid %in% hormone_dat$FSH$IID,
         has_LH = eid %in% hormone_dat$LH$IID)

fsh_dat <- dat_plot %>% filter(has_FSH)
lh_dat <- dat_plot %>% filter(has_LH)

sink(qc_log, append = T)
cat(paste0("# female participants with FSH & age at menopause: ", 
           nrow(fsh_dat), "\n"))
cat(paste0("\t", 
           "# mean (SE) age at menopause: ", 
           mean(fsh_dat$age_at_menopause), " (",  
           sd(fsh_dat$age_at_menopause)/sqrt(nrow(fsh_dat)), ")", "\n"))

cat(paste0("# female participants with LH & age at menopause: ", 
           nrow(lh_dat), "\n"))
cat(paste0("\t", 
           "# mean (SE) age at menopause: ", 
           mean(lh_dat$age_at_menopause), " (",  
           sd(lh_dat$age_at_menopause)/sqrt(nrow(lh_dat)), ")", "\n"))
sink()

# Test for difference in means between groups (with/without FSH & with/without LH) -----

# We want a non-parametric test as age at menopause is not normally distributed
# also try a t-test because it's better powered and the deviation from normality is not huge

no_fsh <- dat_plot %>% filter(!has_FSH)
no_lh <- dat_plot %>% filter(!has_LH)

# fsh
# test_fsh <- wilcox.test(no_fsh$age_at_menopause, fsh_dat$age_at_menopause)
test_fsh <- t.test(no_fsh$age_at_menopause, fsh_dat$age_at_menopause)
# lh
# test_lh <- wilcox.test(no_lh$age_at_menopause, lh_dat$age_at_menopause)
test_lh <- t.test(no_lh$age_at_menopause, lh_dat$age_at_menopause)

sink(qc_log, append = T)
cat(paste0("# P-value for difference between those with & without FSH: ", 
           test_fsh$p.value, "\n"))

cat(paste0("# P-value for difference between those with & without LH: ",
           test_lh$p.value, "\n"))
sink()

# Plot the distribution of age at menopause by FSH/LH status ----

meno_colpalette <- c("#40004B", "#00441B")
names(meno_colpalette) <- c("has_hormone", "no_hormone")

xmin <- min(aam_dat$age_at_menopause)
xmax <- max(aam_dat$age_at_menopause)

no_hr <- dat_plot %>%
  filter(!has_FSH & !has_LH)

p1 <- ggplot(no_hr, aes(x = age_at_menopause)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", alpha = 0.2, binwidth = 1,
                 colour = meno_colpalette["no_hormone"],
                 fill = meno_colpalette["no_hormone"]) +
  geom_vline(xintercept = mean(no_hr$age_at_menopause), linetype = "dashed") +
  scale_x_continuous(limits = c(xmin, xmax)) +
  labs(x = "Age at menopause")
p2 <- ggplot(fsh_dat, aes(x = age_at_menopause)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", alpha = 0.2, binwidth = 1,
                 colour = meno_colpalette["has_hormone"],
                 fill = meno_colpalette["has_hormone"]) +
  geom_vline(xintercept = mean(fsh_dat$age_at_menopause), linetype = "dashed") +
  scale_x_continuous(limits = c(xmin, xmax)) +
  labs(x = "Age at menopause")
p3 <- ggplot(lh_dat, aes(x = age_at_menopause)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", alpha = 0.2, binwidth = 1,
                 colour = meno_colpalette["has_hormone"],
                 fill = meno_colpalette["has_hormone"]) +
  geom_vline(xintercept = mean(lh_dat$age_at_menopause), linetype = "dashed") +
  scale_x_continuous(limits = c(xmin, xmax)) +
  labs(x = "Age at menopause")

age_plot <- ggarrange(p1, p2, p3, nrow = 3)

ggsave(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/plots/aam_distribution.png"), 
       age_plot)

