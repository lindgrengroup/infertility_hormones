# Author: Samvida S. Venkatesh
# Date: 16/01/24

library(tidyverse)
theme_set(theme_bw())

hormones_get <- c("FSH_F", "Oestradiol_F",
                  "Testosterone_F", "Testosterone_M")
genes_get <- read.table("characterise_rare_variants/gene_classification.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

PTHRESH <- 5E-06

gene_res <- lapply(hormones_get, function (hr) {
  res <- read.table(paste0("gene_results/", hr, "_gene.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  
  # unweighted result
  unweighted_res <- read.table(paste0("gene_results/", hr, "_unweighted_gene.txt"),
                               sep = "\t", header = T, stringsAsFactors = F)

  sig_genes <- genes_get %>% filter(strata == hr)
  
  sub_genes <- unique(sig_genes$Region)
  res <- res %>% filter(Region %in% sub_genes) %>%
    filter(Group %in% c("pLoF", "pLoF;damaging_missense_or_protein_altering",
                        "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
                        "Cauchy")) %>%
    select(all_of(c("Region", "Group", "max_MAF",
                    "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                    "BETA_Burden", "SE_Burden"))) 
  
  # Add in unweighted beta
  unweighted_res <- unweighted_res %>% filter(gene_id %in% sub_genes) %>%
    select(all_of(c("gene_id", "Group", "max_MAF",
                    "BETA_Burden_unweighted")))
  
  res <- res %>% left_join(unweighted_res,
                           by = c("Region" = "gene_id",
                                  "Group" = "Group", "max_MAF" = "max_MAF"))
  
  res$hgnc_symbol <- sig_genes$hgnc_symbol[match(res$Region, sig_genes$Region)]
  res$classification <- sig_genes$classification[match(res$Region, sig_genes$Region)]
  return (res)
})
names(gene_res) <- hormones_get

wide_res_tt <- lapply(hormones_get, function (hr) {
  df <- gene_res[[hr]] %>%
    mutate(group_maf = paste0(Group, "_", max_MAF),
           Pvalue = signif(Pvalue, 3),
           Pvalue_Burden = signif(Pvalue_Burden, 3),
           Pvalue_SKAT = signif(Pvalue_SKAT, 3),
           Beta_write_Burden = paste0(signif(BETA_Burden, 3), " (", signif(SE_Burden, 3), ")"),
           Beta_write_Burden_unweighted = signif(BETA_Burden_unweighted, 3))
  
  wide_df <- df %>%
    pivot_wider(id_cols = c("Region", "hgnc_symbol", "classification"),
                names_from = group_maf,
                values_from = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                                "Beta_write_Burden", "Beta_write_Burden_unweighted"))
  for_tt <- wide_df %>% select(any_of(c("classification", "Region", "hgnc_symbol",
                        "Pvalue_Cauchy_NA",
                        "Beta_write_Burden_pLoF_1e-04", "Beta_write_Burden_unweighted_pLoF_1e-04", 
                        "Pvalue_Burden_pLoF_1e-04", "Pvalue_SKAT_pLoF_1e-04",
                        "Beta_write_Burden_pLoF;damaging_missense_or_protein_altering_1e-04", "Beta_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering_1e-04", 
                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering_1e-04", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering_1e-04",
                        "Beta_write_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04", "Beta_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04",
                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04",
                        "Beta_write_Burden_pLoF_0.001", "Beta_write_Burden_unweighted_pLoF_0.001", 
                        "Pvalue_Burden_pLoF_0.001", "Pvalue_SKAT_pLoF_0.001",
                        "Beta_write_Burden_pLoF;damaging_missense_or_protein_altering_0.001", "Beta_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering_0.001", 
                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering_0.001", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering_0.001",
                        "Beta_write_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001", "Beta_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001", 
                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001",
                        "Beta_write_Burden_pLoF_0.01", "Beta_write_Burden_unweighted_pLoF_0.01", 
                        "Pvalue_Burden_pLoF_0.01", "Pvalue_SKAT_pLoF_0.01", 
                        "Beta_write_Burden_pLoF;damaging_missense_or_protein_altering_0.01", "Beta_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering_0.01", 
                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering_0.01", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering_0.01",
                        "Beta_write_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01", "Beta_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01", 
                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01"
  )))
  
  for_tt <- for_tt %>%
    # based on decode feedback only retain cauchy P < 5E-06
    filter(Pvalue_Cauchy_NA < PTHRESH) %>% 
    arrange(classification, Region)
  
  write.table(for_tt, paste0("gene_results/sig_gene_for_manuscript_", hr, ".txt"),
              sep = "\t", row.names = F, quote = F)
  
})

### For infertility ----

infert_get <- c("female_infertility_binary",
                "idiop_infertility_exclusion_binary",
                "male_infertility_binary")
PTHRESH <- 5E-06

gene_res <- lapply(infert_get, function (hr) {
  res <- read.table(paste0("gene_results/", hr, "_gene.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  
  # unweighted result
  unweighted_res <- read.table(paste0("gene_results/", hr, "_unweighted_gene.txt"),
                               sep = "\t", header = T, stringsAsFactors = F)
  
  sig_genes <- res %>%
    mutate(min_p = pmin(Pvalue, Pvalue_Burden, Pvalue_SKAT, na.rm = T)) %>%
    group_by(Region) %>%
    summarise(min_p = min(min_p, na.rm = T)) %>%
    filter(min_p <= PTHRESH)

  res <- res %>% filter(Region %in% unique(sig_genes$Region)) %>%
    filter(Group %in% c("pLoF", "pLoF;damaging_missense_or_protein_altering",
                        "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
                        "Cauchy")) %>%
    select(all_of(c("Region", "Group", "max_MAF",
                    "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                    "BETA_Burden", "SE_Burden"))) 
  
  # Add in unweighted beta
  unweighted_res <- unweighted_res %>% 
    select(all_of(c("gene_id", "Group", "max_MAF",
                    "BETA_Burden_unweighted")))
  
  res <- res %>% left_join(unweighted_res,
                           by = c("Region" = "gene_id",
                                  "Group" = "Group", "max_MAF" = "max_MAF"))
  
  return (res)
})
names(gene_res) <- infert_get

wide_res_tt <- lapply(infert_get, function (hr) {
  df <- gene_res[[hr]] %>%
    mutate(group_maf = paste0(Group, "_", max_MAF),
           Pvalue = signif(Pvalue, 3),
           Pvalue_Burden = signif(Pvalue_Burden, 3),
           Pvalue_SKAT = signif(Pvalue_SKAT, 3),
           OR_Burden = exp(BETA_Burden), 
           LCI_Burden = exp(BETA_Burden - 1.96*SE_Burden), UCI_Burden = exp(BETA_Burden + 1.96*SE_Burden),
           OR_write_Burden = paste0(signif(OR_Burden, 3), " (", signif(LCI_Burden, 3), " - ", signif(UCI_Burden, 3), ")"),
           OR_write_Burden_unweighted = signif(exp(BETA_Burden_unweighted), 3))
  
  wide_df <- df %>%
    pivot_wider(id_cols = Region,
                names_from = group_maf,
                values_from = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                                "OR_write_Burden", "OR_write_Burden_unweighted"))
  for_tt <- wide_df %>% select(any_of(c("Region", 
                                        "Pvalue_Cauchy_NA",
                                        "OR_write_Burden_pLoF_1e-04", "OR_write_Burden_unweighted_pLoF_1e-04", 
                                        "Pvalue_Burden_pLoF_1e-04", "Pvalue_SKAT_pLoF_1e-04",
                                        "OR_write_Burden_pLoF;damaging_missense_or_protein_altering_1e-04", "OR_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering_1e-04", 
                                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering_1e-04", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering_1e-04",
                                        "OR_write_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04", "OR_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04", 
                                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_1e-04",
                                        "OR_write_Burden_pLoF_0.001", "OR_write_Burden_unweighted_pLoF_0.001",
                                        "Pvalue_Burden_pLoF_0.001", "Pvalue_SKAT_pLoF_0.001",
                                        "OR_write_Burden_pLoF;damaging_missense_or_protein_altering_0.001", "OR_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering_0.001", 
                                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering_0.001", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering_0.001",
                                        "OR_write_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001", "OR_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001", 
                                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.001",
                                        "OR_write_Burden_pLoF_0.01", "OR_write_Burden_unweighted_pLoF_0.01",
                                        "Pvalue_Burden_pLoF_0.01", "Pvalue_SKAT_pLoF_0.01",
                                        "OR_write_Burden_pLoF;damaging_missense_or_protein_altering_0.01", "OR_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering_0.01", 
                                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering_0.01", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering_0.01",
                                        "OR_write_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01", "OR_write_Burden_unweighted_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01", 
                                        "Pvalue_Burden_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01", "Pvalue_SKAT_pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous_0.01"
  )))
  
  for_tt <- for_tt %>%
    # based on decode feedback only retain cauchy P < 5E-06
    filter(Pvalue_Cauchy_NA < PTHRESH)
  
  write.table(for_tt, paste0("gene_results/sig_gene_for_manuscript_", hr, ".txt"),
              sep = "\t", row.names = F, quote = F)
  
})















