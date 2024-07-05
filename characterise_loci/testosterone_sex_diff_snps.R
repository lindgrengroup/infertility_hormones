# Author: Samvida S. Venkatesh
# Date: 05/07/2024

library(tidyverse)
theme_set(theme_bw())

# Read data ----

testo_dat <- read.table("Testosterone_all_sexes_gws_results.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

# Bernabeu et al. sex-differential SNPs across the phenome, LD-pruned
phenome_sdsnps <- read.table("number_sd_snps_phenome.txt",
                             sep = "\t", header = T, stringsAsFactors = F,
                             quote = "")

# Prune to get independent lead SNPs across analyses ----

PTHRESH <- 5E-08

all_snps <- testo_dat %>%
  filter(PVALUE_F <= PTHRESH | PVALUE_M <= PTHRESH) %>%
  # filter(AF_Tested_F >= 0.01 & AF_Tested_F <= 0.99) %>%
  # filter(AF_Tested_M >= 0.01 & AF_Tested_M <= 0.99) %>%
  mutate(PVALUE = pmin(PVALUE_F, PVALUE_M))

all_snps <- split(all_snps, f = all_snps$CHROM)

KB_PRUNE <- 500
WINDOW_SIZE <- KB_PRUNE*1000 # 500kb window

pruneSNPs <- function (df, dist_bp = 500000) {
  # Start with the first SNP in the dataframe
  old_top <- -Inf
  current_top <- df$GENPOS[1]
  # Check for SNP with lowest p-value in this region
  while (current_top != old_top) {
    old_top <- current_top
    current_dat <- df %>% filter(GENPOS >= (old_top - dist_bp) &
                                   GENPOS <= (old_top + dist_bp))
    current_top <- current_dat$GENPOS[which.min(current_dat$PVALUE)]
  }
  return (df[df$GENPOS == current_top, ])
}

# Apply to all genome-wide sig SNPs
pruned_snps_full <- lapply(all_snps, function (cdf) {
  print(unique(cdf$CHROM))
  dat_to_prune <- cdf
  pruned_list <- c()
  dist_prune <- KB_PRUNE*1000
  # Repeat until there are no loci remaining
  while (nrow(dat_to_prune) > 0) {
    res_pruned <- pruneSNPs(dat_to_prune, dist_bp = WINDOW_SIZE)
    pruned_list <- c(pruned_list, res_pruned$ID)
    # Remove all SNPs within vicinity of pruned
    dat_to_prune <- dat_to_prune %>%
      filter(GENPOS <= (res_pruned$GENPOS - WINDOW_SIZE) |
               GENPOS >= (res_pruned$GENPOS + WINDOW_SIZE))
  }
  # Results data
  pruned_dat <- cdf %>%
    filter(ID %in% pruned_list)
  return (pruned_dat)
})
pruned_snps_full <- bind_rows(pruned_snps_full)

# Effect sizes for shared SNPs only ----

plot_df <- pruned_snps_full %>%
  filter((PVALUE_F <= 5E-6 & PVALUE_M <= 5E-8) | (PVALUE_F <= 5E-8 & PVALUE_M <= 5E-6)) %>%
  mutate(LCI_M = BETA_M - 1.96*SE_M, UCI_M = BETA_M + 1.96*SE_M, 
         LCI_F = BETA_F - 1.96*SE_F, UCI_F = BETA_F + 1.96*SE_F) 

# Set consistent X and Y axes
max_axis <- max(c(plot_df$UCI_M, plot_df$UCI_F))
min_axis <- min(c(plot_df$LCI_M, plot_df$LCI_F))

p <- ggplot(plot_df, aes(x = BETA_M, y = BETA_F)) +
  geom_pointrange(aes(xmin = LCI_M, xmax = UCI_M),
                  size = 0.3, fatten = 2) +
  geom_pointrange(aes(ymin = LCI_F, ymax = UCI_F),
                  size = 0.3, fatten = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(min_axis, max_axis)) +
  scale_y_continuous(limits = c(min_axis, max_axis)) +
  #geom_abline(intercept = 0, slope = 1) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        plot.title = element_text(size = 8))

ggsave("gws_snps_both_scatter.png", p,
       height = 5, width = 5, units = "cm")


# Identify sex-diff SNPs (Bernabeu et al. define as heterogeneity P < 1E-08) ----

sd_snps <- pruned_snps_full %>%
  filter(het_pval < 1E-08)

# Count on the autosomes and X-chr

autosome_testo_sd <- length(sd_snps$ID[sd_snps$CHROM != 23])
xchr_testo_sd <- length(sd_snps$ID[sd_snps$CHROM == 23])

# Add to phenome-wide table

phenome_sdsnps <- bind_rows(phenome_sdsnps,
                            data.frame(Trait.code = "venkatesh",
                                       Description = "testosterone",
                                       autosomal_sd = autosome_testo_sd,
                                       xchr_sd = xchr_testo_sd))

# Plot ----

# Autosomes

plot_df <- phenome_sdsnps %>%
  mutate(type = ifelse(Trait.code == "venkatesh", "testosterone", "other")) %>%
  filter(autosomal_sd > 0) %>%
  arrange(autosomal_sd)
pheno_levels <- plot_df$Trait.code
plot_df <- plot_df %>%
  mutate(Trait.code = factor(Trait.code, levels = rev(pheno_levels)))

col_palette <- c("#A4A4A4", "#000000")
names(col_palette) <- c("other", "testosterone")

p <- ggplot(data = plot_df, aes(x = Trait.code, 
                                y = autosomal_sd, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_palette, guide = "none") +
  #scale_colour_manual(values = col_palette, guide = "none") +
  labs(y = "Number of sex-differential lead SNPs") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())
ggsave("autosomal_sdsnps.png", p,
       height = 5, width = 15, units = "cm")

# X-chromosome

plot_df <- phenome_sdsnps %>%
  mutate(type = ifelse(Trait.code == "venkatesh", "testosterone", "other")) %>%
  filter(xchr_sd > 0) %>%
  arrange(xchr_sd)
pheno_levels <- plot_df$Trait.code
plot_df <- plot_df %>%
  mutate(Trait.code = factor(Trait.code, levels = rev(pheno_levels)))

col_palette <- c("#A4A4A4", "#000000")
names(col_palette) <- c("other", "testosterone")

p <- ggplot(data = plot_df, aes(x = Trait.code, 
                                y = xchr_sd, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_palette, guide = "none") +
  #scale_colour_manual(values = col_palette, guide = "none") +
  labs(y = "Number of sex-differential lead SNPs") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())
ggsave("xchr_sdsnps.png", p,
       height = 5, width = 15, units = "cm")

# Cross-sex rG ----

# Compare to Bernabeu et al. traits

phenome_cross_sex_rg <- read.table("cross_sex_rg_phenome_wide.txt",
                                   sep = "\t", header = T, stringsAsFactors = F,
                                   quote = "")

plot_df <- phenome_cross_sex_rg %>%
  mutate(type = ifelse(trait_code == "venkatesh", "testosterone", "other")) %>%
  arrange(rg) %>%
  mutate(uci = rg + 1.96*se,
         lci = rg - 1.96*se)
pheno_levels <- plot_df$trait_code
plot_df <- plot_df %>%
  mutate(trait_code = factor(trait_code, levels = rev(pheno_levels)))

col_palette <- c("#A4A4A4", "#000000")
names(col_palette) <- c("other", "testosterone")

p <- ggplot(data = plot_df, aes(y = trait_code, 
                                x = rg, fill = type)) +
  geom_bar(stat = "identity") +
  #geom_errorbar(aes(ymin = lci, ymax = uci)) +
  scale_fill_manual(values = col_palette, guide = "none") +
  #scale_colour_manual(values = col_palette, guide = "none") +
  labs(y = "Genetic correlation") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank())
ggsave("check_ukb_cross_sex_rg.png", p,
       height = 15, width = 7.5, units = "cm")




