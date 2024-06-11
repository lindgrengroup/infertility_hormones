# Author: Samvida S. Venkatesh
# Date: 17/10/23

library(tidyverse)
library(coloc)

INFERT_TRAITS <- c("female_infertility_analysis1", "female_infertility_analysis3",
            "female_infertility_analysis4", "female_infertility_analysis5")
REPRO_TRAITS <- c("endometriosis", "hmb", "PCOS", "uterine_fibroids")

mainpath <- "/well/lindgren/samvida/hormones_infertility/colocalisation"

# Read data ----

lead_snps <- lapply(INFERT_TRAITS, function (infert) {
  res <- read.table(paste0(mainpath, "/infertility_sentinels/",
                    infert, "_eur_sentinels.txt"),
             sep = "\t", header = F, stringsAsFactors = F)
  colnames(res) <- c("MarkerID", "chr", "pos")
  return (res)
})
lead_snps <- bind_rows(lead_snps) %>% distinct()

# Load GWAS results 

infert_res <- lapply(INFERT_TRAITS, function (infert) {
  per_marker <- lapply(lead_snps$MarkerID, function (snp_id) {
    fname <- paste0(mainpath, "/sentinel_windows/", snp_id, "_",
                    infert, ".txt")
    if (file.exists(fname)) {
      res <- read.table(fname, 
                        sep = "\t", header = T, stringsAsFactors = F)
      res$Allele1 <- toupper(res$Allele1)
      res$Allele2 <- toupper(res$Allele2)
      res$markerMatch <- paste0("chr", res$chr, "_", res$pos)
    } else {
      res <- NULL
    }
    return (res)
  })
  names(per_marker) <- lead_snps$MarkerID
  return (per_marker)
})
names(infert_res) <- INFERT_TRAITS

repro_res <- lapply(REPRO_TRAITS, function (repro_tt) {
  per_marker <- lapply(lead_snps$MarkerID, function (snp_id) {
    fname <- paste0(mainpath, "/sentinel_windows/", snp_id, "_",
                    repro_tt, ".txt")
    if (file.size(fname) > 0) {
      res <- read.table(paste0(mainpath, "/sentinel_windows/", snp_id, "_",
                               repro_tt, ".txt"), 
                        sep = "\t", header = T, stringsAsFactors = F)
      res$markerMatch <- paste0(res$CHR, "_", res$POS)
    } else {
      res <- NULL
    }
    return (res)
  })
  names(per_marker) <- lead_snps$MarkerID
  return (per_marker)
})
names(repro_res) <- REPRO_TRAITS

# Wrangle data ----

# For each reproductive trait and marker, create summary stats with matching variants

full_res <- lapply(INFERT_TRAITS, function (infert) {
  per_repro <- lapply(REPRO_TRAITS, function (repro_tt) {
    per_marker <- lapply(lead_snps$MarkerID, function (snp_id) {
      
      inf_df <- infert_res[[infert]][[snp_id]]
      repro_df <- repro_res[[repro_tt]][[snp_id]]
      
      if (!is.null(inf_df) & !is.null(repro_df)) {
        
        # Wrangle infertility results retaining only common >1% markers
        inf_df <- inf_df %>% filter(maf >= 0.01) %>%
          mutate(varbeta_infert = StdErr^2,
                 N_infert = N_CASES + N_CONTROLS,
                 casefrac_infert = N_CASES/N_infert) %>%
          rename(MarkerName_infert = MarkerName,
                 A1_infert = Allele1, A2_infert = Allele2,
                 MAF_infert = maf, BETA_infert = Effect, 
                 PVAL_infert = P.value) %>%
          select(all_of(c("MarkerName_infert", "A1_infert", "A2_infert",
                          "MAF_infert", "BETA_infert", "varbeta_infert",
                          "PVAL_infert", "casefrac_infert", "N_infert",
                          "markerMatch")))
        
        # Wrangle reproductive trait results retaining only common >1% markers
        repro_df <- repro_df %>% 
          mutate(MAF = ifelse(FREQ1 > 0.5, 1-FREQ1, FREQ1),
                 varbeta_repro = SE^2,
                 N_repro = NCASE + NCONTROL,
                 casefrac_repro = NCASE/N_repro) %>%
          filter(MAF >= 0.01) %>%
          rename(MarkerName_repro = RSID,
                 A1_repro = ALLELE1, A2_repro = ALLELE2,
                 MAF_repro = MAF, BETA_repro = BETA, 
                 PVAL_repro = PVALUE) %>%
          select(all_of(c("MarkerName_repro", "A1_repro", "A2_repro",
                          "MAF_repro", "BETA_repro", "varbeta_repro", 
                          "PVAL_repro", "casefrac_repro", "N_repro", 
                          "markerMatch")))
        
        # Merge 
        full_df <- inner_join(inf_df, repro_df, by = "markerMatch")
      } else {
        full_df <- NULL
      }
      return (full_df)
    })
    names(per_marker) <- lead_snps$MarkerID
    return (per_marker)
  })
  names(per_repro) <- REPRO_TRAITS
  return(per_repro)
})
names(full_res) <- INFERT_TRAITS 

# Functions to apply colocalisation ----

map_hypothesis_desc <- data.frame(hypothesis = c("H0", "H1", "H2", "H3", "H4"),
                                  hypo_desc = c("no_assoc", 
                                                "only_infert", "only_repro", 
                                                "both_diff_causal", "both_same_causal"))

applyColoc <- function (df) {
  d_infert <- list(snp = df$markerMatch, 
                   pvalues = df$PVAL_infert, 
                   beta = df$BETA_infert, 
                   varbeta = df$varbeta_infert, 
                   MAF = df$MAF_infert, 
                   s = df$casefrac_infert,
                   N = df$N_infert, type = "cc")
  
  d_repro <- list(snp = df$markerMatch, 
                   pvalues = df$PVAL_repro, 
                   beta = df$BETA_repro, 
                   varbeta = df$varbeta_repro, 
                   MAF = df$MAF_repro, 
                   s = df$casefrac_repro,
                   N = df$N_repro, type = "cc")
  
  coloc_df <- coloc.abf(dataset1 = d_infert, 
                        dataset2 = d_repro, p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)
  
  # Extract and append results
  res <- data.frame(nsnps = coloc_df$summary[1],
                    PP.H0 = coloc_df$summary[2],
                    PP.H1 = coloc_df$summary[3],
                    PP.H2 = coloc_df$summary[4],
                    PP.H3 = coloc_df$summary[5],
                    PP.H4 = coloc_df$summary[6])
  res$ml_hypothesis <- paste0("H", which.max(res[, 2:6])-1)
  return (res)
}

# Apply to all traits ----

all_coloc_res <- lapply(INFERT_TRAITS, function (infert) {
  per_repro <- lapply(REPRO_TRAITS, function (repro_tt) {
    per_marker <- lapply(lead_snps$MarkerID, function (snp_id) {
      df_for_coloc <- full_res[[infert]][[repro_tt]][[snp_id]]
      if (!is.null(df_for_coloc)) {
        res_df <- applyColoc(df_for_coloc)
        res_df$infertility <- infert
        res_df$repro_trait <- repro_tt
        res_df$lead_variant <- snp_id
      } else {
        res_df <- NULL
      }
      return (res_df)
    })
    per_marker <- bind_rows(per_marker)
    return (per_marker)
  })
  per_repro <- bind_rows(per_repro)
  return (per_repro)
})
all_coloc_res <- bind_rows(all_coloc_res)
all_coloc_res$ml_hypothesis_desc <- map_hypothesis_desc$hypo_desc[match(all_coloc_res$ml_hypothesis,
                                                                        map_hypothesis_desc$hypothesis)]

write.table(all_coloc_res,
            paste0(mainpath, "/repro_infertility_coloc_results.txt"),
            sep = "\t", row.names = F, quote = F)

