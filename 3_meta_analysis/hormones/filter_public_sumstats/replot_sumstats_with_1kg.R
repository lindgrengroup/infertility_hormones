# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(argparse)
library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

# Read in arguments ----

parser <- ArgumentParser()
parser$add_argument("--inputFile", required=TRUE,
                    help = "Path to summary statistics")
parser$add_argument("--logFile", required=TRUE,
                    help = "Path to log file to store # SNPs cleaned")
parser$add_argument("--outputFile", required = TRUE,
                    help = "Output file name for filtered summary statistics")
parser$add_argument("--outPlotDir", required = TRUE,
                    help = "Path to store output Manhattan and QQ-plots")
args <- parser$parse_args()

dir.create(args$outPlotDir)

# Wrangle data for cleaning ----

gwas_dat <- read.table(args$inputFile, header = T, 
                       comment.char = "@", stringsAsFactors = F)

# Prepare columns for cleaning
to_numeric <- c("CHROM", "GENPOS", "A1FREQ", "MAF", 
                "BETA", "SE", "STAT", "N", "PVALUE")

gwas_dat <- gwas_dat %>% as_tibble() %>%
  mutate(across(all_of(to_numeric), as.numeric)) 

sink(args$logFile, append = T)
cat(paste0("# SNPs matched 1000 Genomes alleles: ", 
           nrow(gwas_dat), "\n"))
sink()

# Double-check MAF ----

cleaned <- gwas_dat %>%
  mutate(MAF = ifelse(A1FREQ < 0.5, A1FREQ, 1-A1FREQ)) %>%
  # Remove markers with MAF = 0
  filter(MAF > 0)

sink(args$logFile, append = T)
cat(paste0("\t", 
           "# SNPs removed with MAF = 0: ", 
           nrow(gwas_dat) - nrow(cleaned), "\n"))
sink()

# Report metrics post QC
sink(args$logFile, append = T)
cat(paste0("# SNPs post-QC: ", nrow(cleaned), "\n"))
sink()

# Print file formatted for METAL ----

to_print <- cleaned[, c("ID", "CHROM", "GENPOS", 
                        "ALLELE1", "ALLELE0", "A1FREQ", "MAF", 
                        "BETA", "SE", "STAT", "N", "PVALUE")]
# Make sure integers as printed as integers and not in scientific
options(scipen = 999)
write.table(to_print, 
            args$outputFile,
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
  ggsave(paste0(args$outPlotDir, "/with_1kg_MAF_qq_plots.png"),
         full_plot, units = "in", height = 7, width = 7)
} else {
  full_plot <- ggarrange(plotlist = qq_per_bin, nrow = 2, ncol = 3)
  ggsave(paste0(args$outPlotDir, "/with_1kg_MAF_qq_plots.png"),
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

ggsave(paste0(args$outPlotDir, "/with_1kg_MAF_manhattan.png"),
       width = 14, height = 7, units = "in", manhattan_plot)
