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
                      STRATA, "_gene.txt")
outplot_dir <- paste0("/well/lindgren/samvida/hormones_infertility/exome_seq_results/plots/",
                      STRATA)

dir.create(outplot_dir)

two_col_scheme <- c("#202020", "#D35C79")
names(two_col_scheme) <- c("pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
                           "pLoF;damaging_missense_or_protein_altering")

four_col_scheme <- c("#202020", "#53747A" , "#C7B241", "#D35C79")
names(four_col_scheme) <- c("synonymous",
                             "other_missense_or_protein_altering",
                             "damaging_missense_or_protein_altering",
                             "pLoF")

# Read in gene annotations ----

# gene_annots <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/ukb_wes_450k.qced.brava.v1.worst_csq_by_gene_canonical.txt",
#                           sep = "\t", header = T, stringsAsFactors = F)
# 
# gene_annots <- gene_annots %>%
#   mutate(CHROM = gsub(":.*", "", MarkerID),
#          gene_symbol = gsub("/.*", "", gene_symbol),
#          gene_id = gsub("/.*", "", gene_id))
# gene_annots$CHROM <- gsub("^chr", "", gene_annots$CHROM)
# 
# gene_annots <- gene_annots %>%
#   distinct(CHROM, gene_id, gene_symbol)
# 
# gene_annots <- gene_annots[, c("CHROM", "gene_id", "gene_symbol")]
# write.table(gene_annots,
#             "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/gene_canonical.txt",
#             sep = "\t", row.names = F, quote = F)

gene_annots <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/gene_canonical.txt",
                          sep = "\t", header = T, stringsAsFactors = F)

# Wrangle data for cleaning ----

gene_res <- read.table(infile_name, sep = "\t", header = T, 
                       comment.char = "@", stringsAsFactors = F)

# Merge in gene annotations
gene_res <- gene_res %>%
  left_join(gene_annots, by = c("Region" = "gene_id"))

# Convert chrom X to 23
gene_res$CHROM[gene_res$CHROM == "X"] <- 23

# Prepare columns for cleaning
to_numeric <- c("CHROM", "max_MAF", 
                "Pvalue", "Pvalue_SKAT", "Pvalue_Burden", 
                "BETA_Burden", "SE_Burden",
                "MAC", "Number_rare", "Number_ultra_rare",
                "CHROM")

gene_res <- gene_res %>% as_tibble() %>%
  mutate(across(all_of(to_numeric), as.numeric)) 

# Helper functions ----

# Subset data based on MAF and annotation sets
getSubDat <- function (max_maf_cut = 0.001, 
                       test_type = "SKAT-O",
                       annot_type = "combined") {
  # Subset results based on MAF
  sub_res <- gene_res %>% filter(max_MAF == max_maf_cut)
  # Subset results based on annotation category
  if (annot_type == "combined") {
    sub_res <- sub_res %>% 
      filter(Group %in% c("pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
                          "pLoF;damaging_missense_or_protein_altering"))
  } else if (annot_type == "per_category") {
    sub_res <- sub_res %>% 
      filter(Group %in% c("synonymous",
                          "other_missense_or_protein_altering",
                          "damaging_missense_or_protein_altering",
                          "pLoF"))
  }
  
  cols_get <- c("CHROM", "Region", "gene_symbol", 
                "Group", "max_MAF")
  
  if (test_type == "SKAT-O") {
    cols_get <- c(cols_get, "Pvalue")
  } else if (test_type == "SKAT") {
    cols_get <- c(cols_get, "Pvalue_SKAT")
  } else if (test_type == "Burden") {
    cols_get <- c(cols_get, "Pvalue_Burden")
  }
  
  sub_res <- sub_res %>%
    select(all_of(cols_get))
  colnames(sub_res) <- c("CHROM", "Region", "gene_symbol", 
                         "Group", "max_MAF", "Pvalue")
  
  sub_res <- sub_res %>%
    filter(!is.na(CHROM) & !is.na(Pvalue) & Pvalue > 0)
  return (sub_res)
}

# QQ plots and lambdaGC in each MAF bin, coloured by consequence grouping ----

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

getQQ <- function (max_maf_cut = 0.001,
                   maf_bin_name = "MAF < 0.1%", 
                   test_type = "SKAT-O", annot_type = "combined") {
  
  dat_for_plot <- getSubDat(max_maf_cut = max_maf_cut,
                            test_type = test_type, 
                            annot_type = annot_type)
  
  if (annot_type == "combined") { 
    colpal_use <- two_col_scheme 
  } else if (annot_type == "per_category") { 
    colpal_use <- four_col_scheme 
  }
  
  # Get lambdaGC
  chisq_use <- qchisq(1 - dat_for_plot$Pvalue, 1)
  lambdaGC <- median(chisq_use) / qchisq(0.5, 1)
  
  # Split by annotation category to get obs. vs expected pval
  dat_for_plot <- split(dat_for_plot, f = dat_for_plot$Group)
  obs_exp_list <- lapply(dat_for_plot, function (sub_df) {
    qq_df <- getNegLogObsExp(sub_df$Pvalue)
    qq_df$Group <- sub_df$Group[1]
    return (qq_df)
  })
  for_qq_plot <- bind_rows(obs_exp_list)
  
  # For QQ-plot confidence interval ribbon, just subset to one group
  # otherwise we are plotting too many ribbons
  gp_sub <- for_qq_plot$Group[which.max(for_qq_plot$exp)]
  plot_res <- ggplot(for_qq_plot, 
                     aes(x = exp, y = obs)) +
    geom_ribbon(data = for_qq_plot %>% filter(Group == gp_sub),
                aes(ymin = lci, ymax = uci), 
                fill = "grey", alpha = 0.5) +
    geom_point(aes(fill = Group, colour = Group)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = colpal_use) +
    scale_fill_manual(values = colpal_use) +
    labs(title = paste0(maf_bin_name, ", test: ", test_type,
                        ", lambdaGC = ", round(lambdaGC,3)),
         x = "Expected -log10(P)", y = "Observed -log10(P)")
  return (plot_res)
}

# Plot all the QQ plots for the same MAF together
maf_cuts <- sort(unique(gene_res$max_MAF))
names(maf_cuts) <- paste0("MAF <", maf_cuts*100, "%")

qq_plots <- lapply(1:length(maf_cuts), function (i) {
  comb_list <- lapply(c("SKAT-O", "SKAT", "Burden"), function (tty) {
    return (getQQ(max_maf_cut = maf_cuts[i], 
                  maf_bin_name = names(maf_cuts)[i], 
                  test_type = tty, annot_type = "combined"))
  })
  sep_list <- lapply(c("SKAT-O", "SKAT", "Burden"), function (tty) {
    return (getQQ(max_maf_cut = maf_cuts[i], 
                  maf_bin_name = names(maf_cuts)[i], 
                  test_type = tty, annot_type = "per_category"))
  })
  plist <- append(sep_list, comb_list)
  full_plot <- ggarrange(plotlist = plist,
                         ncol = 3, nrow = 2,
                         common.legend = T)
  ggsave(paste0(outplot_dir, "/gene_based_qq_max_maf_", 
                maf_cuts[i], ".png"),
         full_plot, units = "in", height = 7, width = 14)
})

# Cauchy QQ-plot that combines all tests together
cauchy_res <- gene_res %>% filter(Group == "Cauchy")
# Get lambdaGC
chisq_use <- qchisq(1 - cauchy_res$Pvalue, 1)
lambdaGC_cauchy <- median(chisq_use) / qchisq(0.5, 1)

cauchy_qq <- getNegLogObsExp(cauchy_res$Pvalue)

cauchy_qq <- ggplot(cauchy_qq, 
                   aes(x = exp, y = obs)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), 
              fill = "grey", alpha = 0.5) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = paste0("Test: Cauchy, lambdaGC = ", 
                      round(lambdaGC_cauchy, 3)),
       x = "Expected -log10(P)", y = "Observed -log10(P)")

ggsave(paste0(outplot_dir, "/gene_based_qq_cauchy.png"),
       cauchy_qq, units = "in", height = 7, width = 7)

# Manhattan plots ----

gws_two_col_scheme <- c("#A4A4A4", "#A6E8F5", "#202020", "#D35C79")
names(gws_two_col_scheme) <- c("odd_nonsig", "even_nonsig", 
                               "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
                               "pLoF;damaging_missense_or_protein_altering")

gws_four_col_scheme <- c("#A4A4A4", "#A6E8F5",
                         "#202020", "#53747A" , "#C7B241", "#D35C79")
names(gws_four_col_scheme) <- c("odd_nonsig", "even_nonsig", 
                                "synonymous",
                                "other_missense_or_protein_altering",
                                "damaging_missense_or_protein_altering",
                                "pLoF")


getManhattan <- function (max_maf_cut = 0.001,
                          maf_bin_name = "MAF < 0.1%", 
                          test_type = "SKAT-O", annot_type = "combined") {
  # Get data subset
  dat_for_plot <- getSubDat(max_maf_cut = max_maf_cut,
                            test_type = test_type, 
                            annot_type = annot_type)
  
  if (annot_type == "combined") { 
    colpal_use <- gws_two_col_scheme 
  } else if (annot_type == "per_category") { 
    colpal_use <- gws_four_col_scheme 
  }
  
  # Assign a gene position to each gene on the chromosome
  dat_for_plot <- dat_for_plot %>% 
    filter(CHROM %in% 1:23) %>%
    # get chromosome length
    group_by(CHROM) %>% 
    mutate(GENPOS = row_number())
  
  # Add other columns for plotting
  dat_for_plot <- dat_for_plot %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(GENPOS)) %>%
    # get chromosome position
    mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
    # add to original results dataset
    left_join(dat_for_plot, ., by = c("CHROM" = "CHROM")) %>%
    # add cumulative position of each SNP
    arrange(CHROM, GENPOS) %>% 
    mutate(BP_pos = GENPOS + tot) %>%
    # Add highlight and annotation information
    # Add highlight and annotation information
    mutate(status = ifelse(CHROM %% 2 == 0 & Pvalue >= 5e-6, "even_nonsig",
                           ifelse(CHROM %% 2 != 0 & Pvalue >= 5e-6, "odd_nonsig",
                                  ifelse(Pvalue < 5e-6, Group, NA))))
  
  # Axis should just show chromosome number
  axisdf <- dat_for_plot %>% group_by(CHROM) %>% 
    summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)
  
  # Plot
  plot_res <- ggplot(dat_for_plot, aes(x = BP_pos, y = -log10(Pvalue)),
                     fill = status, colour = status) +
    geom_point(data = dat_for_plot %>% filter(status %in% c("even_nonsig", "odd_nonsig")), 
               aes(fill = status, colour = status), 
               shape = 19, size = 1) +
    geom_point(data = dat_for_plot %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
               aes(fill = status, colour = status), 
               shape = 19, size = 1.5) +
    geom_text_repel(data = dat_for_plot %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
                    aes(label = gene_symbol, colour = status), 
                    show.legend = F) + 
    geom_hline(yintercept = -log10(5e-6), linetype = "dashed") +
    scale_colour_manual(values = colpal_use, guide = "none") +
    scale_fill_manual(values = colpal_use, guide = "none") +
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$centre) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = paste0(maf_bin_name, ", test: ", test_type),
         x = "CHROM", y = "-log10(PVAL)") +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  
  return (plot_res)
}

# Split Manhattan plots by MAF
manh_plots <- lapply(1:length(maf_cuts), function (i) {
  comb_list <- lapply(c("SKAT-O", "SKAT", "Burden"), function (tty) {
    p <- getManhattan(max_maf_cut = maf_cuts[i], 
                      maf_bin_name = names(maf_cuts)[i], 
                      test_type = tty, annot_type = "combined")
    ggsave(paste0(outplot_dir, "/gene_based_manhattan_max_maf_", 
                  maf_cuts[i], "_test_", tty, "_annot_combined.png"),
           p, units = "in", height = 7, width = 14)
  })
  sep_list <- lapply(c("SKAT-O", "SKAT", "Burden"), function (tty) {
    p <- getManhattan(max_maf_cut = maf_cuts[i], 
                      maf_bin_name = names(maf_cuts)[i], 
                      test_type = tty, annot_type = "per_category")
    ggsave(paste0(outplot_dir, "/gene_based_manhattan_max_maf_", 
                  maf_cuts[i], "_test_", tty, "_annot_per_category.png"),
           p, units = "in", height = 7, width = 15)
  })
})

# Cauchy Manhattan that combines all tests together
gws_cauchy_col_scheme <- c("#A4A4A4", "#A6E8F5", "#D35C79")
names(gws_cauchy_col_scheme) <- c("odd_nonsig", "even_nonsig", 
                               "signif")
cauchy_res <- gene_res %>% filter(Group == "Cauchy")

# Assign a gene position to each gene on the chromosome
dat_for_plot <- cauchy_res %>% 
  filter(CHROM %in% 1:23) %>%
  # get chromosome length
  group_by(CHROM) %>% 
  mutate(GENPOS = row_number())

# Add other columns for plotting
dat_for_plot <- dat_for_plot %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(GENPOS)) %>%
  # get chromosome position
  mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(dat_for_plot, ., by = c("CHROM" = "CHROM")) %>%
  # add cumulative position of each SNP
  arrange(CHROM, GENPOS) %>% 
  mutate(BP_pos = GENPOS + tot) %>%
  # Add highlight and annotation information
  # Add highlight and annotation information
  mutate(status = ifelse(CHROM %% 2 == 0 & Pvalue >= 5e-6, "even_nonsig",
                         ifelse(CHROM %% 2 != 0 & Pvalue >= 5e-6, "odd_nonsig",
                                ifelse(Pvalue < 5e-6, "signif", NA))))

# Axis should just show chromosome number
axisdf <- dat_for_plot %>% group_by(CHROM) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

# Plot
plot_res <- ggplot(dat_for_plot, aes(x = BP_pos, y = -log10(Pvalue)),
                   fill = status, colour = status) +
  geom_point(data = dat_for_plot %>% filter(status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status), 
             shape = 19, size = 1) +
  geom_point(data = dat_for_plot %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = status, colour = status), 
             shape = 19, size = 1.5) +
  geom_text_repel(data = dat_for_plot %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
                  aes(label = gene_symbol, colour = status), 
                  show.legend = F) + 
  geom_hline(yintercept = -log10(5e-6), linetype = "dashed") +
  scale_colour_manual(values = gws_cauchy_col_scheme, guide = "none") +
  scale_fill_manual(values = gws_cauchy_col_scheme, guide = "none") +
  scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = "Cauchy test",
       x = "CHROM", y = "-log10(PVAL)") +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(outplot_dir, "/gene_based_manhattan_cauchy.png"),
       plot_res, units = "in", height = 7, width = 15)
