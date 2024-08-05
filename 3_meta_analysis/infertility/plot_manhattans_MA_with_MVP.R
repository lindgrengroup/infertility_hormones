# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(argparse)
library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

# Read in arguments ----

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Trait to plot")
args <- parser$parse_args()
STRATA <- args$strata

input_path_eur <- "/well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/results"
input_path_all <- "/well/lindgren/samvida/hormones_infertility/infertility_meta_mvp/filtered"
plot_dir <- "/well/lindgren/samvida/hormones_infertility/infertility_meta_plots"

# Wrangle data for cleaning ----

# EUR results
eur_res <- read.table(paste0(input_path_eur, "/", STRATA, "_eur1.out"),
                      sep = "\t", header = T, stringsAsFactors = F,
                      comment.char = "@")

# All-anc results
all_anc_res <- read.table(paste0(input_path_all, "/", STRATA, "_all_filtered.txt"),
                          sep = "\t", header = T, stringsAsFactors = F,
                          comment.char = "@")

# Standardise columns to plot
to_numeric <- c("Freq1", "FreqSE", "Effect", "StdErr", "P.value",
                "HetISq", "HetChiSq", "HetDf", "HetPVal",
                "GENPOS", "MAF", "BETA", "SE", "PVALUE")

eur_res <- eur_res %>% as_tibble() %>%
  mutate(across(any_of(to_numeric), as.numeric))

all_anc_res <- all_anc_res  %>% as_tibble() %>%
  mutate(across(any_of(to_numeric), as.numeric))

# Extract chromosome and position from MarkerName column 

cleaned_eur <- eur_res %>%
  mutate(chrpos = gsub(":[[:alpha:]].*[[:alpha:]]", "", MarkerName),
         CHROM = gsub(":.*", "", chrpos),
         GENPOS = as.numeric(gsub(".*:", "", chrpos)))
# Replace chromosome X with 23
cleaned_eur$CHROM[which(cleaned_eur$CHROM == "X")] <- 23
cleaned_eur$CHROM <- as.numeric(cleaned_eur$CHROM)
# Add MAF
cleaned_eur <- cleaned_eur %>%
  mutate(MAF = ifelse(Freq1 < 0.5, Freq1, 1-Freq1)) %>%
  rename(ID = "MarkerName", SE = "StdErr")

# Cleaning functions ----

# Only keep autosomes and X chromosome
chrom_filter <- function (dat) {
  res <- dat %>% filter(CHROM %in% c(1:23))
  return (res)
}

# Implausibly large standard error (> 10)
extreme_effect <- function (dat) {
  res <- dat %>% filter(SE < 10)
  return (res)
}

# Remove all instances of markers with duplicate entries
duplicate_snps <- function (dat) {
  res <- dat %>% 
    filter(!(duplicated(ID) | duplicated(ID, fromLast = T)))
  return (res)
}

# Apply cleaning functions ----

cleaned_eur <- chrom_filter(cleaned_eur)
cleaned_eur <- extreme_effect(cleaned_eur)
cleaned_eur <- duplicate_snps(cleaned_eur)

cleaned_all <- chrom_filter(cleaned_all)
cleaned_all <- extreme_effect(cleaned_all)
cleaned_all <- duplicate_snps(cleaned_all)

# Rename columns to standardise
cleaned_eur <- cleaned_eur %>%
  select(any_of(c("MarkerName", "CHROM", "GENPOS", "MAF",
                  "Allele1", "Allele2", "Freq1", "FreqSE", 
                  "Effect", "StdErr", "P.value", "Direction", "HetPVal"))) %>%
  rename(any_of(c(ID = "MarkerName",  
                  BETA = "Effect", SE = "StdErr", PVALUE = "P.value")))

cleaned_all <- cleaned_all %>%
  select(any_of(c("MarkerName", "CHROM", "GENPOS", "MAF",
                  "Allele1", "Allele2", "Freq1", "FreqSE", 
                  "Effect", "StdErr", "P.value", "Direction", "HetPVal"))) %>%
  rename(any_of(c(ID = "MarkerName",  
                  BETA = "Effect", SE = "StdErr", PVALUE = "P.value")))

# Get list of markers that are ONLY significant in all-ancestry analyses
# but not in EUR
all_anc_markers <- cleaned_all$ID[cleaned_all$PVALUE <= 5E-08]
eur_markers <- cleaned_eur$ID[cleaned_eur$PVALUE <= 5E-08]
only_all_markers <- all_anc_markers[!all_anc_markers %in% eur_markers]

# Manhattan plots ----

col_palette <- c("#A4A4A4", "#A6E8F5", "#D35C79", "#009593")
names(col_palette) <- c("odd_nonsig", "even_nonsig", 
                        "signif", "signif_low_MAF")

sub_gwas <- cleaned_all

sub_gwas <- cleaned_all %>% 
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
                                       ifelse(PVALUE < 5e-8 & MAF < 0.01, "signif_low_MAF", NA)))),
         anc_marker = ifelse(ID %in% only_all_markers, "only_all_anc", "eur"))

# Axis should just show chromosome number
axisdf <- sub_gwas %>% group_by(CHROM) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

manhattan_plot <- ggplot(sub_gwas, aes(x = BP_pos, y = -log10(PVALUE)),
                         fill = status, colour = status) +
  geom_point(data = sub_gwas %>% filter(status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status), shape = 19, size = 0.7) +
  geom_point(data = sub_gwas %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status, shape = anc_marker), size = 1) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_colour_manual(values = col_palette, guide = "none") +
  scale_fill_manual(values = col_palette, guide = "none") +
  scale_shape_manual(values = c(only_all_anc = 17, eur = 19), guide = "none") +
  scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8, family = "Arial"))

ggsave(paste0(plot_dir, "/", STRATA, "_manhattan.png"),
       width = 15, height = 4, units = "cm", manhattan_plot)
