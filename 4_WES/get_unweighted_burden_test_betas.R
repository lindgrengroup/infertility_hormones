# Author: Samvida S. Venkatesh
# Date: 16/07/24
# Adapted from: Nik Baya

library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("--phenotype", required=TRUE,
                    help = "Phenotype")
args <- parser$parse_args()

PHENO <- args$phenotype

mainpath <- "/well/lindgren/samvida/hormones_infertility/exome_seq_results/results_2310"

# Read data ----

variant_res <- read.table(paste0(mainpath, "/cat_results/", PHENO,
                                 "_variant.txt"),
                          comment.char = "@", 
                          sep = "\t", header = T, stringsAsFactors = F)

# Variant annotations

variant_annots <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/ukb_wes_450k.july.qced.brava.v6.worst_csq_by_gene_canonical.tsv",
                             sep = "\t", header = T, stringsAsFactors = F)
colnames(variant_annots) <- c("MarkerID", "gene_id", "gene_symbol",
                              "most_severe_consequence", "annotation")

# Wrangle gene ID and variant consequence into variant file ----

# remove duplicates in variant-based results
variant_res <- variant_res %>% distinct(MarkerID, .keep_all = T)

# Merge in annotations, create MAF and align Tstat to minor allele
non_UR_res <- variant_res %>%
  filter(CHR != "UR") %>%
  left_join(variant_annots, by = "MarkerID",
            relationship = "one-to-many") %>%
  mutate(MAF = ifelse(AF_Allele2 > 0.5, 1-AF_Allele2, AF_Allele2),
         Tstat_minor = ifelse(AF_Allele2 > 0.5, -Tstat, Tstat))

# For variants with "UR", the gene id, consequence, MAF is in MarkerID
UR_vars <- variant_res %>%
  filter(CHR == "UR") %>%
  separate(MarkerID, into = c("gene_id", "annotation", "MAF"),
           sep = ":") %>%
  mutate(Tstat_minor = Tstat, MAF = as.numeric(MAF))

variant_res <- bind_rows(non_UR_res, UR_vars)

# Get unweighted gene burden by grouping all the variants for each consequence group and maxMAF combination ----

GROUPS <- c("pLoF", 
            "pLoF;damaging_missense_or_protein_altering",
            "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous")

MAF_WINDOWS <- c(0.01, 0.001, 0.0001)

unweighted_res <- lapply(GROUPS, function (gp) {
  gp_list <- unlist(strsplit(gp, ";"))
  per_maf_bin <- lapply(MAF_WINDOWS, function (max_maf) {
    sub_vars <- variant_res %>%
      filter(annotation %in% gp_list & MAF <= max_maf)
    
    # Summarise across each gene
    grouped_vars <- sub_vars %>% group_by(gene_id) %>%
      summarise(Tstat_sum = sum(Tstat_minor, na.rm = T),
                var_sum = sum(var, na.rm = T)) %>%
      ungroup() %>%
      mutate(BETA_Burden_unweighted = Tstat_sum/var_sum)
    
    # Add back info on group and max_MAF
    grouped_vars$Group <- gp
    grouped_vars$max_MAF <- max_maf
    
    return (grouped_vars)
  })
  per_maf_bin <- bind_rows(per_maf_bin)
  return (per_maf_bin)
})
unweighted_res <- bind_rows(unweighted_res)

write.table(unweighted_res, paste0(mainpath, "/unweighted_betas/", PHENO, "_unweighted_gene.txt"),
              sep = "\t", row.names = F, quote = F)
  
