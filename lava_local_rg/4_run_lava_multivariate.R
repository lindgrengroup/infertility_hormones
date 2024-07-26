# Author: Samvida S. Venkatesh
# Date: 18/07/2024
# Adapted from: https://github.com/josefin-werme/LAVA

library(tidyverse)
library(LAVA)

mainpath <- "/well/lindgren/samvida/hormones_infertility/lava_local"

INFERT <- c("female_infertility_analysis1", "female_infertility_analysis3",
            "female_infertility_analysis5")
HORMONES <- c("fsh_f", "testosterone_f", "williams_tsh", "verdiesen_amh")
REPRO <- c("endometriosis", "hmb", "pcos", "uterine_fibroids", "huber_twinning")
OBESITY <- c("body_fat_percentage", "body_mass_index_bmi_",
             "comparative_body_size_at_age_10", 
             "hip_circumference", "waist_circumference",
             "weight_change_compared_with_1_year_ago",
             "weight", "whole_body_fat_mass")

# Read data ----

# summary statistics
input <- process.input(input.info.file = paste0(mainpath, "/input_sumstats.txt"),
                       sample.overlap.file = paste0(mainpath, "/cross_trait_sample_overlap_from_ldsc.txt"),
                       ref.prefix = "/well/lindgren/samvida/Resources/1000Genomes/for_lava/g1000_eur", # 1000 genomes on GrCh37
                       phenos = c(INFERT, HORMONES, REPRO, OBESITY)) 

# locus files
loci <- read.loci(paste0(mainpath, "/locus_files/LAVA_s2500_m25_f1_w200.loci"))

# Loci with multiple local rG 
multi_loc <- read.table(paste0(mainpath, "/loci_for_multivariate_analysis.txt"),
                        sep = "\t", header = T, stringsAsFactors = F)
multi_loc <- multi_loc %>%
  mutate(chrpos = paste0("chr", chrom, ":", start, "-", stop))
NLOC <- nrow(multi_loc)

# Bivariate results 
fnames <- list.files(paste0(mainpath, "/bivariate_results"), pattern = "*_bivar.txt")
full_fnames <- paste0(mainpath, "/bivariate_results/", fnames)
fnames <- fnames[sapply(full_fnames, file.size) > 1] # get rid of empty files
lava_bivar <- lapply(fnames, function (fn) {
  res <- read.table(paste0(mainpath, "/bivariate_results/", fn), 
                    header = T, stringsAsFactors = F)
  return (res)
})
lava_bivar <- bind_rows(lava_bivar)

# Process bivariate results to select the subset of phenotypes we want to test at each locus ----

PTHRESH <- 2.59E-03 # FDR significance as pre-determined

lava_bivar <- lava_bivar %>% 
  filter(chrpos %in% multi_loc$chrpos & 
           p < PTHRESH &
           phen1 %in% c(INFERT, HORMONES, REPRO, OBESITY) &
           phen2 %in% c(INFERT, HORMONES, REPRO, OBESITY))

multi_loc_phenos <- lapply(1:NLOC, function (li) {
  biv_res <- lava_bivar[lava_bivar$chrpos == multi_loc$chrpos[li], ]
  # At each locus, try all three infertility targets
  per_target_phenos <- lapply(INFERT, function (tgt) {
    biv_tgt_res <- biv_res[biv_res$phen1 == tgt | biv_res$phen2 == tgt, ]
    res_row <- NULL
    if (nrow(biv_tgt_res) > 0) {
      assoc_phenos <- unique(c(biv_tgt_res$phen1, biv_tgt_res$phen2))
      test_phenos <- assoc_phenos[!assoc_phenos %in% INFERT]
      # Add the target and other phenos to dataframe
      if (length(test_phenos) > 1) {
        res_row <- data.frame(chrpos = multi_loc$chrpos[li],
                              target_pheno = tgt,
                              test_phenos = paste0(test_phenos, collapse = ","))
      }
    }
    return (res_row)
  })
  per_target_phenos <- bind_rows(per_target_phenos)
  return (per_target_phenos)
})
multi_loc_phenos <- bind_rows(multi_loc_phenos) %>%
  mutate(chrom = gsub(":.*", "", chrpos),
         chrom = as.numeric(gsub("chr", "", chrom)),
         start = gsub("chr.*:", "", chrpos),
         start = gsub("-.*", "", start),
         stop = gsub(".*-", "", chrpos)) 

write.table(multi_loc_phenos, paste0(mainpath, "/phenos_for_multivariate_analysis.txt"),
            sep = "\t", row.names = F, quote = F)

# Run multivariate regression at each locus ----

multivar_res_all <- lapply(1:nrow(multi_loc_phenos), function (li) {
  res_row <- data.frame(chrpos = multi_loc_phenos$chrpos[li])
  loc_to_process <- loci[which(loci$CHR == multi_loc_phenos$chrom[li] & 
                                 loci$START == multi_loc_phenos$start[li] & 
                                 loci$STOP == multi_loc_phenos$stop[li]), ]
  test_phenos <- unlist(strsplit(multi_loc_phenos$test_phenos[li], ","))
  
  locus <- process.locus(loc_to_process, input, phenos = c(multi_loc_phenos$target_pheno[li],
                                                           test_phenos))
  
  # Some loci do not have sensible covariances, results out of bounds, etc. which we want to ignore
  tryCatch(
    {
      # Run multiple regression (don't run bivariate as we already have those)
      mres <- run.multireg(locus, target = multi_loc_phenos$target_pheno[li],
                           only.full.model = T)
      res_row <- as.data.frame(mres)
      res_row$chrpos <- multi_loc_phenos$chrpos[li]
    },
    error = function(e) {
      res_row$fail_reason <- e$message
      return (res_row)
    }
  )
  return (res_row)
})
multivar_res_all <- bind_rows(multivar_res_all)

write.table(multivar_res_all, paste0(mainpath, "/multivariate_lava_results.txt"),
            sep = "\t", row.names = F, quote = F)
