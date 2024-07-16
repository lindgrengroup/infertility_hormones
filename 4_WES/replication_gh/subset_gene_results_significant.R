# Author: Samvida S. Venkatesh
# Date: 16/01/24

library(tidyverse)
theme_set(theme_bw())

hormones_get <- c("FSH.FEMALE", "oestradiol.FEMALE",
                  "testosterone.FEMALE", "testosterone.MALE")
hormones_og <- c("FSH_F", "Oestradiol_F",
                 "Testosterone_F", "Testosterone_M")
names(hormones_og) <- hormones_get

genes_get <- read.table("characterise_rare_variants/gene_classification.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

PTHRESH <- 5E-06

gene_res <- lapply(hormones_get, function (hr) {
  res <- read.table(gzfile(paste0("replication_GH/GNH.44k_exomes.20240703.value_", hr, ".SAS.SAIGE.gene.txt.gz")),
                    sep = "\t", header = T, stringsAsFactors = F)
  
  sig_genes <- genes_get %>% filter(strata == hormones_og[hr])
  
  sub_genes <- unique(sig_genes$Region)
  res <- res %>% filter(Region %in% sub_genes) %>%
    filter(Group %in% c("pLoF", "pLoF;damaging_missense_or_protein_altering",
                        "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
                        "Cauchy")) %>%
    select(all_of(c("Region", "Group", "max_MAF",
                    "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                    "BETA_Burden", "SE_Burden"))) 
  
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
           Beta_write_Burden = paste0(signif(BETA_Burden, 3), " (", signif(SE_Burden, 3), ")"))
  
  wide_df <- df %>%
    pivot_wider(id_cols = c("Region", "hgnc_symbol", "classification"),
                names_from = group_maf,
                values_from = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                                "Beta_write_Burden"))
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
    arrange(classification, Region)
  
  write.table(for_tt, paste0("replication_GH/sig_gene_for_manuscript_", hr, ".txt"),
              sep = "\t", row.names = F, quote = F)
  
})

### For infertility ----

infert_get <- c("female_infertility4.FEMALE",
                  "male_infertility.MALE")
PTHRESH <- 5E-06

GENES_KEEP <- c("ENSG00000196155", "ENSG00000196946")

gene_res <- lapply(infert_get, function (hr) {
  res <- read.table(gzfile(paste0("replication_GH/GNH.44k_exomes.20240703.", hr, ".SAS.SAIGE.gene.txt.gz")),
                    sep = "\t", header = T, stringsAsFactors = F)

  res <- res %>% filter(Region %in% GENES_KEEP) %>%
    filter(Group %in% c("pLoF", "pLoF;damaging_missense_or_protein_altering",
                        "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
                        "Cauchy")) %>%
    select(all_of(c("Region", "Group", "max_MAF",
                    "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                    "BETA_Burden", "SE_Burden"))) 
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
           OR_write_Burden = paste0(signif(OR_Burden, 3), " (", signif(LCI_Burden, 3), " - ", signif(UCI_Burden, 3), ")"))
  
  wide_df <- df %>%
    pivot_wider(id_cols = Region,
                names_from = group_maf,
                values_from = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                                "OR_write_Burden"))
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
  write.table(for_tt, paste0("replication_GH/sig_gene_for_manuscript_", hr, ".txt"),
              sep = "\t", row.names = F, quote = F)
  
})
