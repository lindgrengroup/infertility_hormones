# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(tidyverse)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

# Read files ----

STRATA <- c("FSH_F", 
            "LH_F", 
            "Oestradiol_F", 
            "Testosterone_F", "Testosterone_M")

infile_dir <- "/well/lindgren/samvida/hormones_infertility/combine_wes_gwas/sig_hits"
outplot_dir <- "/well/lindgren/samvida/hormones_infertility/combine_wes_gwas/plots"

# Common variant GWAS
gwas_sig_hits <- lapply(STRATA, function (st) {
  fname <- paste0(infile_dir, "/", st, "_gwas_sig_hits.txt")
  if (file.exists(fname)) {
    res <- read.table(fname,
                      sep = "\t", header = T, stringsAsFactors = F)
    if (nrow(res) > 0) {
      res <- res %>%
        mutate(strata = st, 
               CHROM = as.character(CHROM), 
               GENPOS = as.integer(GENPOS),
               Allele1 = as.character(Allele1),
               Allele2 = as.character(Allele2)) %>%
        filter(CHROM %in% c(1:23, "X"))
    } else res <- NULL
  } else res <- NULL
  return (res)
})
gwas_sig_hits <- bind_rows(gwas_sig_hits)
gwas_sig_hits$type <- "common_variant"

# Rare variant WES
wes_variant_sig_hits <- lapply(STRATA, function (st) {
  fname <- paste0(infile_dir, "/", st, "_wes_variant_sig_hits.txt")
  if (file.exists(fname)) {
    res <- read.table(fname,
                      sep = "\t", header = T, stringsAsFactors = F)
    if (nrow(res) > 0) {
      res <- res %>%
        mutate(strata = st, CHR = as.character(CHR),
               POS = as.integer(POS),
               Allele1 = as.character(Allele1),
               Allele2 = as.character(Allele2)) %>%
        filter(CHR %in% c(1:23, "X"))
    } else res <- NULL
  } else res <- NULL
  return (res)
})
wes_variant_sig_hits <- bind_rows(wes_variant_sig_hits)
wes_variant_sig_hits$type <- "rare_variant"

# Gene-based tests
wes_gene_sig_hits <- lapply(STRATA, function (st) {
  fname <- paste0(infile_dir, "/", st, "_wes_gene_sig_hits.txt")
  if (file.exists(fname)) {
    res <- read.table(fname,
                      sep = "\t", header = T, stringsAsFactors = F)
    if (nrow(res) > 0) {
      res <- res %>%
        mutate(strata = st)
    } else res <- NULL
  } else res <- NULL
  return (res)
})
wes_gene_sig_hits <- bind_rows(wes_gene_sig_hits)
wes_gene_sig_hits$type <- "gene_based"

# Read in variant classification ----

variant_classification <- read.table("/well/lindgren/samvida/hormones_infertility/conditional_analysis/classified_all_lead_snps_all_strata.txt",
                                     sep = "\t", header = T, stringsAsFactors = F)
variant_classification <- variant_classification %>%
  mutate(CHROM = as.character(CHROM))

gene_classification <- read.table("/well/lindgren/samvida/hormones_infertility/characterise_rare_variants/all_strata_gene_classification.txt",
                                  sep = "\t", header = T, stringsAsFactors = F)
gene_classification <- gene_classification %>% distinct()

# Read in gene annotations ----

# Common variant GWAS 
common_variant_annots <- read.table(paste0(infile_dir, "/all_hormones_sig_gwas_annots.txt"),
                                    sep = "\t", header = T, stringsAsFactors = F,
                                    fill = NA)
colnames(common_variant_annots) <- c("CHR", "POS", "CHRPOS", 
                                     "gene_id.fpred_max_label", "fpred_max_score", "hgnc_symbol")
# Retain only one annotation per variant to make plots easier to interpret (max fpred score)
common_variant_annots <- common_variant_annots %>%
  filter(!is.na(fpred_max_score)) %>%
  group_by(CHRPOS) %>%
  slice_max(fpred_max_score)
# Separate gene id from variant consequence
common_variant_annots <- common_variant_annots %>%
  mutate(gene_id = gsub(":.*", "", gene_id.fpred_max_label),
         fpred_max_label = gsub(".*:", "", gene_id.fpred_max_label)) %>%
  select(all_of(c("CHRPOS", "gene_id", "fpred_max_label", "fpred_max_score")))

# Add HGNC symbol to gene id 
source("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/scripts/post_hoc/utils.R")
mapping <- get_mapping_ensembl_to_hgnc()
common_variant_annots$hgnc_symbol <- mapping[common_variant_annots$gene_id]
common_variant_annots <- common_variant_annots %>% distinct() %>%
  group_by(CHRPOS) %>%
  summarise(gene_id = str_c(gene_id, collapse = "; "),
            fpred_max_label = str_c(fpred_max_label, collapse = "; "),
            hgnc_symbol = str_c(hgnc_symbol, collapse = "; "),
            fpred_max_score = str_c(fpred_max_score, collapse = "; "))

# Common variant annots: CHRPOS, gene_id, hgnc_symbol, fpred_max_label, fpred_max_score

# Rare variant WES
rare_variant_annots <- read.table(paste0(infile_dir, "/all_hormones_sig_wes_variant_annots.txt"),
                                  sep = "\t", header = T, stringsAsFactors = F)
colnames(rare_variant_annots) <- c("MarkerID", "gene_id", "gene_symbol",
                                   "most_severe_consequence", "annotation")
# Keep the consequence of most impact
consq_priority <- c("pLoF", "damaging_missense",
                    "other_missense", "synonymous", "non_coding")
rare_variant_annots <- rare_variant_annots %>%
  mutate(annotation = factor(annotation, levels = consq_priority)) %>%
  group_by(MarkerID) %>% arrange(annotation, .by_group = T) %>%
  slice(1)

# Gene-based tests
gene_based_annots <- read.table(paste0(infile_dir, "/all_hormones_sig_wes_gene_annots.txt"),
                                sep = "\t", header = T, stringsAsFactors = F)
# Gene annots: CHROM, gene_id, gene_symbol

# Wrangle data for EFFECT SIZE VS MAF PLOTS ----

# Merge in annotations and synthesise column names to be similar across types of data
# COLUMN NAMES: MarkerID, CHROM, GENPOS, MAF, BETA, SE, PVALUE,
# gene_id, hgnc_symbol, consequence_category
# Don't allow multiple consequences per hit (stratify and take the most consequential)

GENE_COLOURS <- c("#C7B241", "#63A369", "#009593", "#CD875D", "#D35C79",
                           "#202020")
                           
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("novel_any_hormone", 
                                 "novel_this_hormone", "reported")

gwas_sig_hits <- left_join(gwas_sig_hits, common_variant_annots,
                           by = c("ID" = "CHRPOS"))
gwas_sig_beta_maf <- gwas_sig_hits %>% 
  mutate(MarkerID = paste0(ID, ":", Allele1, ":", Allele2),
         fpred_max_label = paste0("v2g_", gsub("_variant", "", fpred_max_label)),
         strata = gsub("_EUR", "", strata)) %>%
  select(all_of(c("MarkerID", "CHROM", "GENPOS", "MAF", "BETA", "SE", "PVALUE",
                  "gene_id", "hgnc_symbol", "fpred_max_label",
                  "strata", "type"))) %>%
  rename(annotation = fpred_max_label) 
gwas_sig_beta_maf <- left_join(gwas_sig_beta_maf, 
                               variant_classification[, c("CHROM", "GENPOS", "classification", "strata")],
                               by = c("CHROM", "GENPOS", "strata"))

wes_variant_sig_hits <- left_join(wes_variant_sig_hits, rare_variant_annots,
                                  by = "MarkerID")
wes_variant_beta_maf <- wes_variant_sig_hits %>%
  mutate(MAF = ifelse(AF_Allele2 > 0.5, 1-AF_Allele2, AF_Allele2)) %>%
  select(all_of(c("MarkerID", "CHR", "POS", "MAF", "BETA", "SE", "p.value",
                  "gene_id", "gene_symbol", "annotation",
                  "strata", "type"))) %>%
  rename(CHROM = CHR, GENPOS = POS, PVALUE = p.value,
         hgnc_symbol = gene_symbol) 

wes_gene_sig_hits <- left_join(wes_gene_sig_hits, gene_based_annots,
                               by = c("Region" = "gene_id"))
wes_gene_sig_beta_maf <- wes_gene_sig_hits %>%
  filter(!is.na(BETA_Burden)) %>%
  mutate(Group = paste0("burden_", Group)) %>%
  select(all_of(c("CHROM", "max_MAF", "BETA_Burden", "SE_Burden", "Pvalue_Burden",
                  "Region", "gene_symbol", "Group",
                  "strata", "type"))) %>%
  rename(MAF = max_MAF, BETA = BETA_Burden, SE = SE_Burden, PVALUE = Pvalue_Burden,
         gene_id = Region, hgnc_symbol = gene_symbol, annotation = Group)
wes_gene_sig_beta_maf <- left_join(wes_gene_sig_beta_maf, 
                                   gene_classification[, c("Region", "strata", "classification")],
                                   by = c("gene_id" = "Region", 
                                          "strata" = "strata"))

# Combine all dfs to plot
for_beta_maf_plots <- bind_rows(bind_rows(gwas_sig_beta_maf, wes_variant_beta_maf),
                                wes_gene_sig_beta_maf) 
for_beta_maf_plots$classification[is.na(for_beta_maf_plots$classification)] <- "reported"

## Function to plot BETA vs MAF ----

getBetaMAFPlotbyGene <- function (df, plot_title) {
  df <- df %>%
    mutate(beta_plot = abs(BETA),
           shape_dirn = ifelse(BETA < 0, "down", "up"))
  
  
  genes_to_colour <- df %>% arrange(desc(beta_plot))
  genes_to_colour <- unique(genes_to_colour$gene_id)
  genes_to_colour <- genes_to_colour[!is.na(genes_to_colour)]
  if (length(genes_to_colour) >= 5) {
    genes_to_colour <- genes_to_colour[1:5]
    names(GENE_COLOURS) <- c(genes_to_colour, "other")
  } else {
    names(GENE_COLOURS) <- c(genes_to_colour, 
                             rep("none", times = 5 - length(genes_to_colour)), 
                             "other")
  }
  df <- df %>%
    mutate(shape_colour = ifelse(gene_id %in% genes_to_colour &!is.na(gene_id), 
                                 gene_id, "other"))
  # Label (at most) one hit for each gene
  label_df <- df %>% filter(gene_id %in% genes_to_colour) %>%
    group_by(hgnc_symbol) %>%
    arrange(desc(beta_plot), .by_group = T) %>%
    slice(1) %>%
    arrange(desc(beta_plot)) %>%
    mutate(hgnc_symbol = str_trunc(hgnc_symbol, 10, "right"))
  colpal_use <- GENE_COLOURS
  
  res_plot <- ggplot(df, aes(x = MAF, y = beta_plot,
                             fill = shape_colour, colour = shape_colour)) +
    geom_point(data = df %>% filter(shape_colour == "other"),
               aes(shape = shape_dirn)) +
    geom_point(data = df %>% filter(shape_colour != "other"),
               aes(shape = shape_dirn)) +
    geom_text_repel(data = label_df, 
                    aes(label = hgnc_symbol)) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    scale_fill_manual(values = colpal_use) +
    scale_colour_manual(values = colpal_use) +
    scale_shape_manual(values = c("up" = 24, "down" = 25)) +
    scale_size(range = c(0, 4)) +
    labs(x = "Minor allele frequency", y = "Effect size",
         title = plot_title) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  return (res_plot)
}

getBetaMAFPlotbyNovelty <- function (df, plot_title) {
  df <- df %>%
    mutate(beta_plot = abs(BETA),
           shape_dirn = ifelse(BETA < 0, "down", "up"))
  
  colpal_use <- custom_three_diverge
  # Label (at most) one hit for each gene
  label_df <- df %>% filter(classification != "reported") %>%
    group_by(hgnc_symbol) %>%
    arrange(desc(beta_plot), .by_group = T) %>%
    slice(1) %>%
    arrange(desc(beta_plot)) %>%
    mutate(hgnc_symbol = str_trunc(hgnc_symbol, 10, "right"))
  
  res_plot <- ggplot(df, aes(x = MAF, y = beta_plot,
                             fill = classification, colour = classification)) +
    geom_point(data = df %>% filter(type == "common_variant" & MAF > 0.01)) +
    geom_point(data = df %>% filter(MAF <= 0.01),
               aes(shape = shape_dirn)) +
    geom_text_repel(data = label_df, 
                    aes(label = hgnc_symbol)) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    scale_fill_manual(values = colpal_use) +
    scale_colour_manual(values = colpal_use) +
    scale_shape_manual(values = c("up" = 24, "down" = 25)) +
    scale_size(range = c(0, 4)) +
    labs(x = "Minor allele frequency", y = "Effect size",
         title = plot_title) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  return (res_plot)
}


## Plot separately within each strata and separately for variants and genes ----

lapply(STRATA, function (st) {
  print(st)
  # Common and rare variants
  df_plot_variants <- for_beta_maf_plots %>%
    filter(strata == st & type != "gene_based")
  
  if (nrow(df_plot_variants) > 0) {
    ggsave(paste0(outplot_dir, "/variants_common_rare_", st, "_coloured_by_gene.png"),
           getBetaMAFPlotbyGene(df_plot_variants,
                          plot_title = paste0(st, ": variants")), 
           units = "in", height = 7, width = 7)
  }
  
  df_plot_burden <- for_beta_maf_plots %>%
    filter(strata == st & type != "rare_variant")
  
  if (nrow(df_plot_burden) > 0) {
    ggsave(paste0(outplot_dir, "/gene_burden_common_variants_", st, "_coloured_by_gene.png"),
           getBetaMAFPlotbyGene(df_plot_burden,
                          plot_title = paste0(st, ": gene-burden & common variants")), 
           units = "in", height = 7, width = 7)
    
    ggsave(paste0(outplot_dir, "/gene_burden_common_variants_", st, "_coloured_by_novelty.png"),
           getBetaMAFPlotbyNovelty(df_plot_burden,
                          plot_title = paste0(st, ": gene-burden & common variants")), 
           units = "in", height = 7, width = 7)
  }
})
