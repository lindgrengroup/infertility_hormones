# Author: Samvida S. Venkatesh
# Date: 19/04/2023

library(tidyverse)

mainpath <- "/well/lindgren/samvida/hormones_infertility/characterise_rare_variants"

HORMONES <- c("FSH", "LH", "Oestradiol", "Testosterone")
SEX_STRATA <- c("F", "M", "sex_comb")

# Read GWAS Catalog mapped genes ----

gwascat_mapped_genes <- read.table("/well/lindgren/samvida/Resources/GWASCatalog/gwascat_hormone_associations.txt",
                                  sep = "\t", comment.char = "|",
                                  header = T, stringsAsFactors = F)
gwascat_mapped_genes <- gwascat_mapped_genes %>%
  select(all_of(c("MAPPED_TRAIT", "MAPPED_GENE"))) %>%
  separate_rows(MAPPED_GENE, sep = "[^[:alnum:]]") %>%
  filter(!MAPPED_GENE %in% c("", "NA")) %>%
  distinct()

gwascat_traits <- lapply(HORMONES, function (hr) {
  traits_check <- read.table(paste0("/well/lindgren/samvida/Resources/GWASCatalog/", hr,
                                    "_traits.txt"),
                             sep = "\t", header = F, stringsAsFactors = F)
  return (traits_check$V1)
})
names(gwascat_traits) <- HORMONES

# Significant variants ----

rare_variant_annots <- read.table("/well/lindgren/samvida/hormones_infertility/combine_wes_gwas/sig_hits/all_hormones_sig_wes_variant_annots.txt",
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

sig_res <- lapply(HORMONES, function (hr) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    fname <- paste0("/well/lindgren/samvida/hormones_infertility/combine_wes_gwas/sig_hits/", 
                    hr, "_", sx, "_wes_variant_sig_hits.txt")
    if (file.exists(fname)) {
      res <- read.table(fname, sep = "\t", header = T, stringsAsFactors = F)
      res <- res %>%
        mutate(across(all_of(c("MarkerID", "CHR", "Allele1", "Allele2")), as.character),
               across(all_of(c("POS", "AC_Allele2", "AF_Allele2", "MissingRate",
                               "BETA", "SE", "Tstat", "var", "p.value", "N")), as.numeric)) %>%
        filter(CHR %in% c(1:23, "X"))
      res <- res %>% left_join(rare_variant_annots, by = "MarkerID")
    } else res <- NULL
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(sig_res) <- HORMONES

# Variant pruning ----

# Only retain the lowest P-value rare variant annotated to each gene

pruned_res <- lapply(HORMONES, function (hr) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- sig_res[[hr]][[sx]]
    res <- NULL
    if (!is.null(df)) {
      res <- df %>% filter(AF_Allele2 <= 0.01 | AF_Allele2 >= 0.99) %>%
        separate_rows(gene_symbol, sep = "/") %>% 
        filter(!gene_symbol %in% c("", "NA")) %>%
        filter(!is.na(gene_symbol)) 
      
      res <- res %>% group_by(gene_symbol) %>%
        slice_min(p.value)
      res$hormone <- hr
      res$sex_strata <- sx
      res$gene_symbol <- as.character(res$gene_symbol)
    }
    return (res)
    })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(pruned_res) <- HORMONES

# Classify as reported gene (from GWAS) or novel ----

ANY_HORMONE_GENES <- unique(gwascat_mapped_genes$MAPPED_GENE)
this_hormone_genes <- lapply(HORMONES, function (hr) {
  res_list <- lapply(gwascat_traits[[hr]], function (htr) {
    return (unique(gwascat_mapped_genes$MAPPED_GENE[grepl(htr, 
                                                          gwascat_mapped_genes$MAPPED_TRAIT)]))
  })
  res_list <- unique(unlist(res_list))
  return (res_list)
})
names(this_hormone_genes) <- HORMONES

classif_pruned_res <- lapply(HORMONES, function (hr) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- pruned_res[[hr]][[sx]]
    df$classification <- ifelse(df$gene_symbol %in% this_hormone_genes[[hr]], "reported",
                                        ifelse(df$gene_symbol %in% ANY_HORMONE_GENES, "novel_this_hormone",
                                               "novel_any_hormone"))
    df$classification <- as.character(df$classification)
    return (df)
  })
  res_list <- bind_rows(res_list)
  return (res_list)
})
classif_pruned_res <- bind_rows(classif_pruned_res)

write.table(classif_pruned_res,
            paste0(mainpath, "/all_strata_pruned_annotated_results.txt"),
            sep = "\t", row.names = F, quote = F)

# Plot the number of reported and novel genes in each strata ----

custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("novel_any_hormone", 
                                 "novel_this_hormone", "reported")
plotBars <- function (df) {
  summ_df <- df %>% 
    group_by(strata) %>%
    count(classification)
  
  resplot <- ggplot(summ_df, 
                    aes(x = strata, y = n)) +
    geom_col(aes(fill = classification), width = 0.7) +
    scale_fill_manual(values = custom_three_diverge) +
    labs(x = "Strata", y = "Number of genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  return (resplot)
}

lapply(HORMONES, function (hr) {
  for_plot <- classif_pruned_res %>% filter(hormone == hr) %>%
    mutate(strata = paste0(hormone, "_", sex_strata))
  
  # Barplots
  ggsave(paste0(mainpath, "/plots/", hr, "_number_genes_novel_reported.png"),
         plotBars(for_plot), units = "in", height = 7, width = 7)
})

