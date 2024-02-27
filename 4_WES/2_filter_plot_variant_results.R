# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(argparse)
library(tidyverse)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

# Read in arguments ----

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Strata being assessed")
args <- parser$parse_args()

STRATA <- args$strata

infile_name <- paste0("/well/lindgren/samvida/hormones_infertility/exome_seq_results/results_2310/cat_results/",
                      STRATA, "_variant.txt")
outplot_dir <- paste0("/well/lindgren/samvida/hormones_infertility/exome_seq_results/plots/",
                      STRATA)

dir.create(outplot_dir)

four_col_scheme <- c("#202020", "#202020", "#53747A", "#C7B241", "#D35C79")
names(four_col_scheme) <- c("non_coding", "synonymous",
                           "other_missense", "damaging_missense", 
                            "pLoF")

# Read in variant annotations ----

variant_annots <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/ukb_wes_450k.july.qced.brava.v6.worst_csq_by_gene_canonical.tsv",
                             sep = "\t", header = T, stringsAsFactors = F)
colnames(variant_annots) <- c("MarkerID", "gene_id", "gene_symbol",
                              "most_severe_consequence", "annotation")

# Wrangle data for cleaning ----

variant_res <- read.table(infile_name, sep = "\t", header = T, 
                          comment.char = "@", stringsAsFactors = F)

# Merge in annotations
variant_res <- variant_res %>%
  left_join(variant_annots, by = "MarkerID")

# Turn chrom X to 23
variant_res$CHR[variant_res$CHR == "X"] <- 23

# Prepare columns for cleaning
to_numeric <- c("CHR", "POS", 
                "AC_Allele2", "AF_Allele2",
                "MissingRate", "BETA", "SE", "Tstat",
                "var", "p.value")

variant_res <- variant_res %>% as_tibble() %>%
  filter(CHR %in% 1:23 & !is.na(p.value) & p.value > 0) %>%
  mutate(across(all_of(to_numeric), as.numeric),
         MAF = ifelse(AF_Allele2 <= 0.5, AF_Allele2, 1-AF_Allele2))

min_bb <- ceiling(log10(min(variant_res$MAF)))
if (min_bb <= -3) {
  log_bin_breaks <- min_bb:-3
  bin_breaks <- c(0, 10^log_bin_breaks, 1)
} else {
  bin_breaks <- c(0, 10^min_bb, 1)
}

if (length(bin_breaks) == 3) {
  bin_labels <- c("<0.1%", ">=0.1%")
} else if (length(bin_breaks) == 4) {
  bin_labels <- c("<0.01%", "[0.01% - 0.1%)", ">=0.1%")
} else if (length(bin_breaks) == 5) {
  bin_labels <- c("<0.001%", "[0.001% - 0.01%)", "[0.01% - 0.1%)", ">=0.1%")
} 

variant_res$MAF_bin <- cut(variant_res$MAF, 
                           breaks = bin_breaks, labels = bin_labels,
                           include.lowest = T)

# QQ plots per MAF bin, coloured by consequence category ----

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
  dat_for_plot <- variant_res %>% 
    filter(MAF_bin == mb)
  
  # Get lambdaGC
  chisq_use <- qchisq(1 - dat_for_plot$p.value, 1)
  lambdaGC <- median(chisq_use) / qchisq(0.5, 1)
  
  # Split by annotation category to get obs. vs expected pval
  dat_for_plot <- split(dat_for_plot, f = dat_for_plot$annotation)
  obs_exp_list <- lapply(dat_for_plot, function (sub_df) {
    qq_df <- getNegLogObsExp(sub_df$p.value)
    qq_df$annotation <- sub_df$annotation[1]
    return (qq_df)
  })
  for_qq_plot <- bind_rows(obs_exp_list)
  
  # For QQ-plot confidence interval ribbon, just subset to one consequence
  # otherwise we are plotting too many ribbons
  cat_sub <- for_qq_plot$annotation[which.max(for_qq_plot$exp)]
  plot_res <- ggplot(for_qq_plot, 
                     aes(x = exp, y = obs)) +
    geom_ribbon(data = for_qq_plot %>% filter(annotation == cat_sub),
                aes(ymin = lci, ymax = uci), 
                fill = "grey", alpha = 0.5) +
    geom_point(aes(fill = annotation, colour = annotation)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = four_col_scheme) +
    scale_fill_manual(values = four_col_scheme) +
    labs(title = paste0(mb, ", lambdaGC = ", round(lambdaGC,3)),
         x = "Expected -log10(P)", y = "Observed -log10(P)")
  return (plot_res)
})

full_plot <- ggarrange(plotlist = qq_per_bin,
                       ncol = 2, nrow = 2,
                       common.legend = T)
ggsave(paste0(outplot_dir, "/variant_based_qq_plots.png"),
       full_plot, units = "in", height = 7, width = 7)

# Manhattan plots ----

gws_four_col_scheme <- c("#A4A4A4", "#A6E8F5", "#202020", "#202020", "#53747A", "#C7B241", "#D35C79")
names(gws_four_col_scheme) <- c("odd_nonsig", "even_nonsig", 
                                "non_coding", "synonymous",
                                "other_missense", "damaging_missense", 
                                "pLoF")

manhattan_per_bin <- lapply(1:length(bin_labels), function (mbi) {
  dat_for_plot <- variant_res %>% 
    filter(MAF_bin == bin_labels[mbi] & !is.na(p.value) & p.value > 0) 
  
  # Add other columns for plotting
  dat_for_plot <- dat_for_plot %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POS)) %>%
    # get chromosome position
    mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
    # add to original results dataset
    left_join(dat_for_plot, ., by = c("CHR" = "CHR")) %>%
    # add cumulative position of each SNP
    arrange(CHR, POS) %>% 
    mutate(BP_pos = POS + tot) %>%
    # Add highlight and annotation information
    mutate(status = ifelse(CHR %% 2 == 0 & p.value >= 1e-7, "even_nonsig",
                           ifelse(CHR %% 2 != 0 & p.value >= 1e-7, "odd_nonsig",
                                  ifelse(p.value < 1e-7, annotation, NA))))
  
  # Axis should just show chromosome number
  axisdf <- dat_for_plot %>% group_by(CHR) %>% 
    summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)
  
  # Label significant points with gene, but only label the top of each chimney
  for_label <- dat_for_plot %>%
    filter(p.value < 1e-7) %>%
    mutate(gene_symbol = gsub("/.*", "", gene_symbol)) %>%
    group_by(annotation, gene_symbol) %>%
    slice_min(order_by = p.value) %>%
    select(all_of(c("CHR", "BP_pos", "p.value",
                    "gene_symbol", "annotation", "status")))
  # If there are more than 10 genes to label, pick the top 10 unique genes
  # per chromosome
  nuniq <- length(unique(for_label$gene_symbol))
  if (nuniq > 10) {
    sub_genes <- for_label %>% 
      group_by(gene_symbol) %>% 
      slice_min(order_by = p.value) %>%
      ungroup() %>% group_by(CHR) %>%
      arrange(p.value) %>% slice(1:10)
    
    keep_genes <- unique(sub_genes$gene_symbol)
    for_label <- for_label %>%
      filter(gene_symbol %in% keep_genes)
  }
  
  # Plot
  plot_res <- ggplot(dat_for_plot, aes(x = BP_pos, y = -log10(p.value)),
                     fill = status, colour = status) +
    geom_point(data = dat_for_plot %>% filter(status %in% c("even_nonsig", "odd_nonsig")), 
               aes(fill = status, colour = status), 
               shape = 19, size = 1) +
    geom_point(data = dat_for_plot %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
               aes(fill = status, colour = status), 
               shape = 19, size = 1.5) +
    geom_text_repel(data = for_label, 
                    aes(label = gene_symbol, colour = status), 
                    show.legend = F) + 
    geom_hline(yintercept = -log10(1e-7), linetype = "dashed") +
    scale_colour_manual(values = gws_four_col_scheme, guide = "none") +
    scale_fill_manual(values = gws_four_col_scheme, guide = "none") +
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = paste0("MAF: ", bin_labels[mbi]), x = "CHROM", y = "-log10(PVAL)") +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  
  ggsave(paste0(outplot_dir, 
                "/variant_based_manhattan_max_maf_", bin_breaks[mbi+1], ".png"),
         plot_res, units = "in", height = 7, width = 15)
})

# Run Manhattan for all variants, not separated by MAF bin
dat_for_plot <- variant_res %>% 
  filter(!is.na(p.value) & p.value > 0) 

# Add other columns for plotting
dat_for_plot <- dat_for_plot %>%
  group_by(CHR) %>%
  summarise(chr_len = max(POS)) %>%
  # get chromosome position
  mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(dat_for_plot, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, POS) %>% 
  mutate(BP_pos = POS + tot) %>%
  # Add highlight and annotation information
  mutate(status = ifelse(CHR %% 2 == 0 & p.value >= 1e-7, "even_nonsig",
                         ifelse(CHR %% 2 != 0 & p.value >= 1e-7, "odd_nonsig",
                                ifelse(p.value < 1e-7, annotation, NA))))

# Axis should just show chromosome number
axisdf <- dat_for_plot %>% group_by(CHR) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

# Label significant points with gene, but only label the top of each chimney
for_label <- dat_for_plot %>%
  filter(p.value < 1e-7) %>%
  mutate(gene_symbol = gsub("/.*", "", gene_symbol)) %>%
  group_by(annotation, gene_symbol) %>%
  slice_min(order_by = p.value) %>%
  select(all_of(c("CHR", "BP_pos", "p.value",
                  "gene_symbol", "annotation", "status")))
# If there are more than 10 genes to label, pick the top 10 unique genes
# per chromosome
nuniq <- length(unique(for_label$gene_symbol))
if (nuniq > 10) {
  sub_genes <- for_label %>% 
    group_by(gene_symbol) %>% 
    slice_min(order_by = p.value) %>%
    ungroup() %>% group_by(CHR) %>%
    arrange(p.value) %>% slice(1:10)
  
  keep_genes <- unique(sub_genes$gene_symbol)
  for_label <- for_label %>%
    filter(gene_symbol %in% keep_genes)
}

# Plot
plot_res <- ggplot(dat_for_plot, aes(x = BP_pos, y = -log10(p.value)),
                   fill = status, colour = status) +
  geom_point(data = dat_for_plot %>% filter(status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status), 
             shape = 19, size = 1) +
  geom_point(data = dat_for_plot %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status), 
             shape = 19, size = 1.5) +
  geom_text_repel(data = for_label, 
                  aes(label = gene_symbol, colour = status), 
                  show.legend = F) + 
  geom_hline(yintercept = -log10(1e-7), linetype = "dashed") +
  scale_colour_manual(values = gws_four_col_scheme, guide = "none") +
  scale_fill_manual(values = gws_four_col_scheme, guide = "none") +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "CHROM", y = "-log10(PVAL)") +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(outplot_dir, 
              "/variant_based_manhattan_all_maf.png"),
       plot_res, units = "in", height = 7, width = 15)
