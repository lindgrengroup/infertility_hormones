# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(argparse)
library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

# Read in arguments ----

parser <- ArgumentParser()
parser$add_argument("--index", required=TRUE,
                    help = "Row number of descriptives sheet to run")
args <- parser$parse_args()

main_filepath <- "/well/lindgren/samvida/hormones_infertility/public_sumstats"

# Read study descriptives and GWAS data ----

all_study_desc <- read.table(paste0(main_filepath, "/public_sumstats_descriptives.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
study_desc <- all_study_desc[as.numeric(args$index), ]
study_desc$lifted_filename <- paste0(study_desc$PHENOTYPE, "/lifted/",
                                     gsub(".log", "", gsub(".*filtered/", "", study_desc$logFileName)),
                                     "_hg38.txt")

GWAS_res <- read.table(paste0(main_filepath, "/", study_desc$lifted_filename), 
                         sep = " ", header = T,
                         stringsAsFactors = F, comment.char = "@")

# Wrangle data for cleaning ----

# Prepare columns for cleaning
to_numeric <- c("GENPOS", "A1FREQ", "MAF", "BETA", "SE",
                "STAT", "N", "PVALUE")

GWAS_res <- GWAS_res %>% as_tibble() %>%
  mutate(across(any_of(to_numeric), as.numeric)) 

if (grepl("^chr", GWAS_res$CHROM[1])) {
  GWAS_res$CHROM <- gsub("^chr", "", GWAS_res$CHROM)
}
# Replace X with 23
GWAS_res$CHROM[which(GWAS_res$CHROM == "X")] <- 23
GWAS_res$CHROM <- as.numeric(GWAS_res$CHROM)

# Cleaning functions ----

# Only keep autosomes and X chromosome
chrom_filter <- function (qc_log, dat) {
  res <- dat %>% filter(CHROM %in% c(1:23))
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed not on autosomes or X chromosome: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# INFO score 
info_filter <- function (qc_log, dat) {
  res <- dat %>% filter(INFO >= 0.8)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with low INFO < 0.8: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Only keep A/G/C/T SNPs (no indels)
acgt_filter <- function (qc_log, dat) {
  res <- dat %>% filter(ALLELE1 %in% c("A", "C", "G", "T") &
                          ALLELE0 %in% c("A", "C", "G", "T"))
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# non-ACGT SNPs removed: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Implausibly large standard error (> 10)
extreme_effect <- function (qc_log, dat) {
  res <- dat %>% filter(SE < 10)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with extreme standard error > 10: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Remove all instances of markers with duplicate entries
duplicate_snps <- function (qc_log, dat) {
  res <- dat %>% 
    mutate(dup_check = paste0(CHROM, ":", GENPOS)) %>%
    filter(!(duplicated(dup_check) | duplicated(dup_check, fromLast = T)))
  
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed due to duplicates: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Apply cleaning functions ----

log_file <- paste0(main_filepath, "/", study_desc$logFileName)

sink(log_file, append = T)
cat(paste0("# SNPs in lifted summary statistics file: ", 
           nrow(GWAS_res), "\n"))
sink()

cleaned <- chrom_filter(log_file, GWAS_res)
if ("INFO" %in% colnames(cleaned)) {
  cleaned <- info_filter(log_file, cleaned)
}
cleaned <- acgt_filter(log_file, cleaned)
cleaned <- extreme_effect(log_file, cleaned)
cleaned <- duplicate_snps(log_file, cleaned)

# Report metrics post QC
sink(log_file, append = T)
cat(paste0("# SNPs post-QC: ", nrow(cleaned), "\n"))
sink()

# Print file formatted for METAL ----

# Sort results by chromosome and position and rename SNP to "chrN:pos"
cleaned <- cleaned %>% arrange(CHROM, GENPOS) %>%
  mutate(ID = paste0("chr", CHROM, ":", GENPOS))
if ("LOG10P" %in% colnames(cleaned)) {
  cleaned <- cleaned %>% 
    mutate(PVALUE = 10^-LOG10P)
}

to_print <- cleaned[, c("ID", "CHROM", "GENPOS", 
                        "ALLELE1", "ALLELE0", "A1FREQ", "MAF", 
                        "BETA", "SE", "STAT", "N", "PVALUE")]
# Make sure integers as printed as integers and not in scientific
options(scipen = 999)
write.table(to_print, 
            paste0(main_filepath, "/", study_desc$outputFileName),
            sep = "\t", row.names = F, quote = F)
options(scipen = 0)

# Get MAF bins ----

gwas_dat <- to_print

log_bin_breaks <- ceiling(log10(min(gwas_dat$MAF))):-1
bin_breaks <- c(0, 10^log_bin_breaks, 1)
if (length(bin_breaks) == 3) {
  bin_labels <- c("<10%", ">=10%")
} else if (length(bin_breaks) == 4) {
  bin_labels <- c("<1%", "[1% - 10%)", ">=10%")
} else if (length(bin_breaks) == 5) {
  bin_labels <- c("<0.1%", "[0.1% - 1%)", "[1% - 10%)", ">=10%")
} else {
  # There are very rare variants
  bin_breaks <- c(0, 0.0001, 0.001, 0.01, 0.1, 1)
  bin_labels <- c("<0.01%", "[0.01% - 0.1%)", "[0.1% - 1%)", "[1% - 10%)", ">=10%")
} 
gwas_dat$MAF_bin <- cut(gwas_dat$MAF, 
                        breaks = bin_breaks, labels = bin_labels,
                        include.lowest = T)

# QQ plots and lambdaGC in each MAF bin ----

dir.create(paste0(main_filepath, "/", study_desc$outPlotDir))

getNegLogObsExp <- function (observed_pvals) {
  # Get expected p-values
  obs_pvals <- sort(observed_pvals)
  neglog_obs <- -log10(obs_pvals)
  exp_pvals <- (1:length(obs_pvals) - 0.5)/length(obs_pvals)
  neglog_exp <- -log10(exp_pvals)
  
  # CI based on the beta distribution
  # parameters: a = k, b = n-k where k=seq(1, n)
  k <- 1:length(obs_pvals)
  beta_lci <- qbeta(0.025, k, length(obs_pvals)-k)
  beta_uci <- qbeta(0.975, k, length(obs_pvals)-k)
  
  # Flip lci and uci for -log10
  plot_qq <- data.frame(obs = neglog_obs, exp = neglog_exp,
                        lci = -log10(beta_uci), uci = -log10(beta_lci))
  return (plot_qq)
}

qq_per_bin <- lapply(bin_labels, function (mb) {
  dat_for_plot <- gwas_dat %>% 
    filter(MAF_bin == mb & !is.na(PVALUE) & PVALUE > 0)
  
  # Get lambdaGC
  chisq_use <- qchisq(1 - dat_for_plot$PVALUE, 1)
  lambdaGC <- median(chisq_use) / qchisq(0.5, 1)
  
  qq_df <- getNegLogObsExp(dat_for_plot$PVALUE)
  
  plot_res <- ggplot(qq_df, 
                     aes(x = exp, y = obs)) +
    geom_ribbon(aes(ymin = lci, ymax = uci), 
                fill = "grey", alpha = 0.5) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = paste0(mb, ", lambdaGC = ", round(lambdaGC,3)),
         x = "Expected -log10(P)", y = "Observed -log10(P)")
  return (plot_res)
})

if (length(qq_per_bin) <= 4) {
  full_plot <- ggarrange(plotlist = qq_per_bin, nrow = 2, ncol = 2)
  ggsave(paste0(main_filepath, "/", study_desc$outPlotDir, "/qq_plots.png"),
         full_plot, units = "in", height = 7, width = 7)
} else {
  full_plot <- ggarrange(plotlist = qq_per_bin, nrow = 2, ncol = 3)
  ggsave(paste0(main_filepath, "/", study_desc$outPlotDir, "/qq_plots.png"),
         full_plot, units = "in", height = 7, width = 14)
} 

# Manhattan plots ----

col_palette <- c("#A4A4A4", "#A6E8F5", "#D35C79", "#009593")
names(col_palette) <- c("odd_nonsig", "even_nonsig", 
                        "signif", "signif_low_MAF")

sub_gwas <- gwas_dat
max_n_plot <- max(gwas_dat$N)

sub_gwas <- sub_gwas %>% 
  # get chromosome length
  group_by(CHROM) %>% summarise(chr_len = as.numeric(max(GENPOS))) %>%
  # get chromosome position
  mutate(tot = as.numeric(cumsum(chr_len) - chr_len)) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(sub_gwas, ., by = c("CHROM" = "CHROM")) %>%
  # add cumulative position of each SNP
  arrange(CHROM, GENPOS) %>% mutate(BP_pos = as.numeric(GENPOS + tot)) %>%
  # Add highlight and annotation information
  mutate(status = ifelse(CHROM %% 2 == 0 & PVALUE >= 5e-8, "even_nonsig",
                         ifelse(CHROM %% 2 != 0 & PVALUE >= 5e-8, "odd_nonsig",
                                ifelse(PVALUE < 5e-8 & MAF >= 0.01, "signif", 
                                       ifelse(PVALUE < 5e-8 & MAF < 0.01, "signif_low_MAF", NA)))))

# Axis should just show chromosome number
axisdf <- sub_gwas %>% group_by(CHROM) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

# Plot
manhattan_plot <- ggplot(sub_gwas, aes(x = BP_pos, y = -log10(PVALUE)),
                         fill = status, colour = status) +
  geom_point(data = sub_gwas %>% filter(status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status), shape = 19, size = 1) +
  geom_point(data = sub_gwas %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status), shape = 19, size = 1.5) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_colour_manual(values = col_palette, guide = "none") +
  scale_fill_manual(values = col_palette, guide = "none") +
  scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = paste0("Max N = ", max_n_plot)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(main_filepath, "/", study_desc$outPlotDir, "/manhattan.png"),
       width = 14, height = 7, units = "in", manhattan_plot)
