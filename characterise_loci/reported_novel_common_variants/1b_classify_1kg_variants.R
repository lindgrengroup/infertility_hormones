# Author: Samvida S. Venkatesh
# Date: 22/06/2022

library(tidyverse)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/conditional_analysis"

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone") 
SEX_STRATA <- c("F", "M", "sex_comb")
ANCESTRY_GROUPS <- c("EUR", "all_anc")

# Read original lead SNPs ----

lead_snps <- lapply(HORMONES, function (hr) {
  res_list_sex <- lapply(SEX_STRATA, function (sx) {
    res_list_anc <- lapply(ANCESTRY_GROUPS, function (anc_gp) {
      fname <- paste0(mainpath, "/lead_snps/1kg_lead_snp_sumstats_", hr, "_", sx, "_", anc_gp, 
                      "_with_rsids.txt")
      res <- NULL
      if (file.exists(fname)) {
        if (length(readLines(fname)) > 1) {
          res <- read.table(fname, sep = "\t", header = T, stringsAsFactors = F)
        } 
      }
      return (res)
    })
    names(res_list_anc) <- ANCESTRY_GROUPS
    return (res_list_anc)
  })
  names(res_list_sex) <- SEX_STRATA
  return (res_list_sex)
})
names(lead_snps) <- HORMONES

# Read known novel SNPs ----

# Convert SNPs to "chrN:pos"
confirmed_novel_snps <- lapply(HORMONES, function (hr) {
  res_list_sex <- lapply(SEX_STRATA, function (sx) {
    res_list_anc <- lapply(ANCESTRY_GROUPS, function (anc_gp) {
      fname_specific <- paste0(mainpath, "/novel_snps/", hr, "_", sx, "_", anc_gp,
                               "_novel_specific.txt")
      if (file.exists(fname_specific) & file.size(fname_specific) > 0) {
        hr_specific_novel <- read.table(fname_specific, 
                                        sep = "\t", header = F, stringsAsFactors = F)$V1
        hr_specific_novel <- paste0("chr", gsub("_", ":", hr_specific_novel))
      } else hr_specific_novel <- NA
      
      fname_any <- paste0(mainpath, "/novel_snps/", hr, "_", sx, "_", anc_gp,
                          "_novel_any.txt")
      if (file.exists(fname_any) & file.size(fname_any) > 0) {
        any_novel <- read.table(fname_any, 
                                sep = "\t", header = F, stringsAsFactors = F)$V1
        any_novel <- paste0("chr", gsub("_", ":", any_novel))
      } else any_novel <- NA
      
      return (list(specific = hr_specific_novel,
                   any = any_novel))
    })
    names(res_list_anc) <- ANCESTRY_GROUPS
    return (res_list_anc)
  })
  names(res_list_sex) <- SEX_STRATA
  return (res_list_sex)
})
names(confirmed_novel_snps) <- HORMONES

# Read and classify results from conditional analyses ----

P_CONDTHRESH <- 0.05
isCondIndpt <- function (fname) {
  # Read the conditional analysis result
  if (file.exists(fname) & file.size(fname) > 0) {
    pc_df <- read.table(fname,
                        sep = "\t", header = T, stringsAsFactors = F)
    # The SNP is novel if it's conditional P < 0.05 against all other published SNPs
    res <- all(pc_df$pC < P_CONDTHRESH)
  } else res <- F
  return (res) 
}

poss_novel_snps <- lapply(HORMONES, function (hr) {
  res_list_sex <- lapply(SEX_STRATA, function (sx) {
    res_list_anc <- lapply(ANCESTRY_GROUPS, function (anc_gp) {
      #cat(paste0(hr, "_", sx, "_", anc_gp, "\n"))
      fname_specific <- paste0(mainpath, "/unreported_snps/", hr, "_", sx, "_", anc_gp,
                               "/", hr, "_poss_novel_snp_list.txt")
      if (file.exists(fname_specific) & file.size(fname_specific) > 0) {
        hr_specific_poss_novel <- read.table(fname_specific, 
                                             sep = "\t", header = F, stringsAsFactors = F)$V1
        
        classified_hr_specific <- sapply(hr_specific_poss_novel, function (v) {
          isCondIndpt(paste0(mainpath, "/unreported_snps/", hr, "_", sx, "_", anc_gp,
                             "/", v, "_cond_on_pub_", hr, ".cma.cojo"))
        })
        
        recreate_ids <- strsplit(hr_specific_poss_novel, ":")
        chrpos_ids <- unlist(lapply(recreate_ids, 
                                    function (vec) paste0("chr", vec[1], ":", vec[2])))
        res_specific <- chrpos_ids[classified_hr_specific]
        
        } else res_specific <- NA
      
      fname_any <- paste0(mainpath, "/unreported_snps/", hr, "_", sx, "_", anc_gp,
                          "/any_poss_novel_snp_list.txt")
      if (file.exists(fname_any) & file.size(fname_any) > 0) {
        hr_any_poss_novel <- read.table(fname_any, 
                                        sep = "\t", header = F, stringsAsFactors = F)$V1
        
        classified_hr_any <- sapply(hr_any_poss_novel, function (v) {
          isCondIndpt(paste0(mainpath, "/unreported_snps/", hr, "_", sx, "_", anc_gp,
                             "/", v, "_cond_on_pub_any.cma.cojo"))
        })
        
        recreate_ids <- strsplit(hr_any_poss_novel, ":")
        chrpos_ids <- unlist(lapply(recreate_ids, 
                                    function (vec) paste0("chr", vec[1], ":", vec[2])))
        res_any <- chrpos_ids[classified_hr_any]
      } else res_any <- NA
      
      return (list(specific = res_specific,
                   any = res_any))
    })
    names(res_list_anc) <- ANCESTRY_GROUPS
    return (res_list_anc)
  })
  names(res_list_sex) <- SEX_STRATA
  return (res_list_sex)
})
names(poss_novel_snps) <- HORMONES

# Annotate finemapped SNPs with categorisation ----

CHAR_COLS <- c("ID", "RSID", "Allele1", "Allele2",
               "Direction", "VARID_1KG", "classification",
               "strata", "hormone", "sex_strata", "ancestry")

annot_snps <- lapply(HORMONES, function (hr) {
  res_list_sex <- lapply(SEX_STRATA, function (sx) {
    res_list_anc <- lapply(ANCESTRY_GROUPS, function (anc_gp) {
      log_file <- paste0(mainpath, "/logs/", hr, "_", sx, "_", anc_gp, ".txt")
      df_to_classify <- lead_snps[[hr]][[sx]][[anc_gp]]
      
      if (!is.null(df_to_classify)) {
        novel_snps_specific_hr <- unique(c(confirmed_novel_snps[[hr]][[sx]][[anc_gp]][["specific"]],
                                           poss_novel_snps[[hr]][[sx]][[anc_gp]][["specific"]]))
        
        novel_snps_any_hr <- unique(c(confirmed_novel_snps[[hr]][[sx]][[anc_gp]][["any"]],
                                      poss_novel_snps[[hr]][[sx]][[anc_gp]][["any"]]))
        
        sink(log_file, append = T)
        cat(paste0("\t", "# Novel associations with this hormone: ",
                   sum(!is.na(novel_snps_specific_hr)), 
                   "\n"))
        cat(paste0("\t", "# Novel associations with any hormone: ",
                   sum(!is.na(novel_snps_any_hr)), 
                   "\n"))
        sink()
        
        # Classify SNPs
        df_to_classify$classification <- ifelse(df_to_classify$ID %in% novel_snps_any_hr, 
                                                "novel_any_hormone",
                                                ifelse(df_to_classify$ID %in% novel_snps_specific_hr, 
                                                       "novel_this_hormone", "reported"))
        
        df_to_classify <- df_to_classify %>%
          mutate(strata = paste0(hr, "_", sx, "_", anc_gp),
                 hormone = hr,
                 sex_strata = sx,
                 ancestry = anc_gp) %>%
          mutate(across(all_of(CHAR_COLS), as.character))
      }
      return (df_to_classify)
    })
    res_list_anc <- bind_rows(res_list_anc)
    return (res_list_anc)
  })
  res_list_sex <- bind_rows(res_list_sex)
  return (res_list_sex)
})
annot_snps <- bind_rows(annot_snps)

write.table(annot_snps,
            paste0(mainpath, "/classified_1KG_lead_snps_all_strata.txt"), 
            sep = "\t", row.names = F, quote = F)
