# Author: Samvida S. Venkatesh
# Date: 17/04/2022

library(tidyverse)

mainpath <- "/well/lindgren/samvida/hormones_infertility/conditional_analysis"

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone") 
SEX_STRATA <- c("F", "M", "sex_comb")
ANCESTRY_GROUPS <- c("EUR", "all_anc")

# Make sure columns are in the right types
CHAR_COLS <- c("ID", "RSID", "Allele1", "Allele2", 
               "Direction", "VARID_1KG", "hormone", "sex_strata", "ancestry",
               "classification", "strata")
NUM_COLS <- c("CHROM", "GENPOS", "MAF", "Freq1", "FreqSE",
              "BETA", "SE", "PVALUE", "HetPVal")

# Read lead SNPs ----

lead_snps <- lapply(HORMONES, function (hr) {
  res_list_sex <- lapply(SEX_STRATA, function (sx) {
    res_list_anc <- lapply(ANCESTRY_GROUPS, function (anc_gp) {
      fname <- paste0(mainpath, "/lead_snps/all_lead_snp_sumstats_", hr, "_", sx, "_", anc_gp, 
                      "_with_rsids.txt")
      res <- NULL
      if (file.exists(fname)) {
        if (length(readLines(fname)) > 1) {
          res <- read.table(fname, sep = "\t", header = T,
                            stringsAsFactors = F)
          res <- res %>%
            mutate(across(any_of(CHAR_COLS), as.character)) %>%
            mutate(across(any_of(NUM_COLS), as.numeric)) %>%
            mutate(sex_strata = recode(sex_strata, "FALSE" = "F")) 
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

# Classified 1KG SNPs ----

classified_1kg_snps <- read.table(paste0(mainpath, "/classified_1KG_lead_snps_all_strata.txt"),
                                  sep = "\t", header = T, stringsAsFactors = F)

classified_1kg_snps <- classified_1kg_snps %>%
  mutate(across(any_of(CHAR_COLS), as.character)) %>%
  mutate(across(any_of(NUM_COLS), as.numeric)) %>%
  mutate(sex_strata = recode(sex_strata, "FALSE" = "F")) 

# Reported variant lists ----

reported_variants <- lapply(HORMONES, function (hr) {
  res <- read.table(paste0(mainpath, "/published_snps/", hr, "_reported.txt"),
                    sep = " ", header = F, stringsAsFactors = F)
  colnames(res) <- c("ID", "VARID_1KG", "reported_trait")
  res <- res %>%
    mutate(chromosome = gsub(":.*", "", ID),
           position = gsub(".*:", "", ID))
  res$chromosome <- gsub("chr", "", res$chromosome)
  res$chromosome[res$chromosome == "X"] <- 23
  return (res[, c("ID", "chromosome", "position")])
})
names(reported_variants) <- HORMONES

reported_variants_all_hormones <- read.table(paste0(mainpath, "/published_snps/all_hormones_reported.txt"),
                                             sep = " ", header = F, stringsAsFactors = F)
colnames(reported_variants_all_hormones) <- c("ID", "VARID_1KG", "reported_trait")
reported_variants_all_hormones <- reported_variants_all_hormones %>%
  mutate(chromosome = gsub(":.*", "", ID),
         position = gsub(".*:", "", ID))
reported_variants_all_hormones$chromosome <- gsub("chr", "", reported_variants_all_hormones$chromosome)
reported_variants_all_hormones$chromosome[reported_variants_all_hormones$chromosome == "X"] <- 23
reported_variants_all_hormones <- 
  reported_variants_all_hormones[, c("ID", "chromosome", "position")]

# Functions for sorting ----

WINDOW_SIZE <- 500000 # 500kb window

checkVarinRange <- function (var_to_sort, mid_of_ranges_to_check, 
                             dist_bp = WINDOW_SIZE) {
  var_to_sort <- as.numeric(var_to_sort)
  mid_of_ranges_to_check <- as.numeric(mid_of_ranges_to_check)
  res <- any(sapply(mid_of_ranges_to_check, function (mp) {
    var_to_sort >= (mp-dist_bp) & var_to_sort <= (mp+dist_bp)
  }))
  return (res)
}

# Flag and remove any variants already classified in 1KG
sortVariants <- function (hr, sx, anc_gp) {

  df_to_class <- lead_snps[[hr]][[sx]][[anc_gp]]
  
  if (!is.null(df_to_class)) {
    relevant_vars_1kg <- classified_1kg_snps %>%
      filter(strata == paste0(hr, "_", sx, "_", anc_gp))
    
    # Check by chromosome
    classed_df <- lapply(1:23, function (chrom) {
      # cat(paste0("Testing chromosome: ", chrom, "\n"))
      # Remove 1KG variants
      sub_df_to_class <- df_to_class %>% filter(CHROM == chrom)
      sub_1kg_check <- relevant_vars_1kg %>% filter(CHROM == chrom)
      
      if (nrow(sub_df_to_class) > 0) {
        # Default assume no proxy
        sub_df_to_class$has_1kg_proxy <- F
        has_proxy <- NULL
        if (nrow(sub_1kg_check) > 0) {
          sub_df_to_class <- sub_df_to_class %>%
            rowwise() %>%
            mutate(has_1kg_proxy = checkVarinRange(GENPOS, 
                                                   sub_1kg_check$GENPOS))
          # HAS 1KG PROXY
          has_proxy <- sub_df_to_class %>% 
            filter(has_1kg_proxy) %>% select(-has_1kg_proxy) %>%
            mutate(classification = "has_1kg_proxy")
        }
        sub_df_to_class <- sub_df_to_class %>% 
          filter(!has_1kg_proxy)
      } else has_proxy <- NULL
      
      # Now check whether there is a nearby published variant for *any* hormone
      sub_reported_variants_all_hormones <- reported_variants_all_hormones %>%
        filter(chromosome == chrom)
      
      if (nrow(sub_df_to_class) > 0) {
        if (nrow(sub_reported_variants_all_hormones) > 0) {
          sub_df_to_class <- sub_df_to_class %>%
            rowwise() %>%
            mutate(published_any = checkVarinRange(GENPOS, 
                                                   sub_reported_variants_all_hormones$position),
                   classification = ifelse(published_any, "novel_this_hormone", 
                                           "novel_any_hormone"))
          # NOVEL SNPs
          novel_snps <- sub_df_to_class %>% 
            filter(classification == "novel") %>%
            select(-any_of(c("published_specific",
                             "published_any",
                             "has_1kg_proxy")))
          
        } else if (nrow(sub_df_to_class) > 0) {
          sub_df_to_class$classification <- "novel_any_hormone"
          # NOVEL SNPs
          novel_snps <- sub_df_to_class %>% 
            mutate(classification = "novel_any_hormone") 
        } else novel_snps <- NULL
        sub_df_to_class <- sub_df_to_class %>% 
          filter(classification != "novel_any_hormone")
      } else novel_snps <- NULL
      
      # Now check whether there is a nearby published variant for this specific hormone
      sub_reported_variants <- reported_variants[[hr]] %>%
        filter(chromosome == chrom)
      
      if (nrow(sub_df_to_class) > 0) {
        if (nrow(sub_reported_variants) > 0) {
          sub_df_to_class <- sub_df_to_class %>%
            rowwise() %>%
            mutate(published_specific = checkVarinRange(GENPOS, 
                                                        sub_reported_variants$position),
                   classification = ifelse(published_specific, "reported", 
                                           classification))
        }
        published_res <- sub_df_to_class %>%
          select(-any_of(c("published_specific",
                           "published_any",
                           "has_1kg_proxy")))
      } else published_res <- NULL
      
      # Bind all results
      chr_res <- bind_rows(has_proxy, novel_snps,
                           published_res)
      
      return (chr_res)
    })
    classed_df <- bind_rows(classed_df) %>%
      mutate(across(any_of(CHAR_COLS), as.character)) %>%
      mutate(across(any_of(NUM_COLS), as.numeric))%>%
      mutate(sex_strata = recode(sex_strata, "FALSE" = "F")) 
  } else classed_df <- NULL
  return (classed_df)
}

# Apply to all strata ----

classified_all_snps <- lapply(HORMONES, function (hr) {
  per_sex_strata <- lapply(SEX_STRATA, function (sx) {
    per_anc_group <- lapply(ANCESTRY_GROUPS, function (anc) {
      
      snps_from_1kg <- classified_1kg_snps %>%
        filter(strata == paste0(hr, "_", sx, "_", anc))
      
      # Filter out SNPs with proxies and bind rows with 1KG classification
      non_1kg_snp_classification <- sortVariants(hr, sx, anc)
      if (!is.null(non_1kg_snp_classification)) {
        non_1kg_snp_classification <- non_1kg_snp_classification %>%
          filter(classification != "has_1kg_proxy")
      }

      # Bind rows
      snps_classified <- bind_rows(non_1kg_snp_classification,
                                   snps_from_1kg)
      
      if (!is.null(snps_classified)) {
        snps_classified <- snps_classified %>% 
          arrange(CHROM, GENPOS) %>%
          mutate(across(any_of(CHAR_COLS), as.character)) %>%
          mutate(across(any_of(NUM_COLS), as.numeric)) %>%
          mutate(sex_strata = recode(sex_strata, "FALSE" = "F")) 
        
        # Log results 
        log_file <- paste0(mainpath, "/logs/all_snp_classification_log.txt")
        sink(log_file, append = T)
        cat(paste0("STRATA: ", hr, "_", sx, "_", anc, "\n",
                   "\t", "# Total independent loci: ",
                   nrow(snps_classified),
                   "\n",
                   "\t", "# Novel loci: ",
                   length(which(snps_classified$classification == "novel_any_hormone")),
                   "\n",
                   "\t", "# Previously reported to associate with ", hr, ": ",
                   length(which(snps_classified$classification == "reported")),
                   "\n",
                   "\t", "# Previously reported to associate with SOME hormone: ",
                   length(which(snps_classified$classification == "novel_this_hormone")),
                   "\n"))
        sink()
      }
      return (snps_classified)
    })
    per_anc_group <- bind_rows(per_anc_group)
    return (per_anc_group)
  })
  per_sex_strata <- bind_rows(per_sex_strata)
  return (per_sex_strata)
  })
classified_all_snps <- bind_rows(classified_all_snps) %>%
  mutate(across(any_of(CHAR_COLS), as.character)) %>%
  mutate(across(any_of(NUM_COLS), as.numeric)) %>%
  mutate(sex_strata = recode(sex_strata, "FALSE" = "F"),
         strata = paste0(hormone, "_", sex_strata, "_", ancestry)) 

write.table(classified_all_snps,
            paste0(mainpath, "/classified_all_lead_snps_all_strata.txt"), 
            sep = "\t", row.names = F, quote = F)
