# Author: Samvida S. Venkatesh
# Date: 08/07/2024
# Adapted from: https://github.com/josefin-werme/LAVA/blob/main/vignettes/sample_overlap.Rmd

library(tidyverse)

# Read in tables of LDSC cross-trait rG

# between hormones and infertility
hi_cor <- read.table("/well/lindgren/samvida/hormones_infertility/lava_local/hormones_infertility.rg",
                     header = T, stringsAsFactors = F)
hi_cor <- hi_cor[ ,c("p1","p2","gcov_int")]

# between hormones or infertility and repro diseases
repro_cor <- read.table("/well/lindgren/samvida/hormones_infertility/lava_local/repro_traits.rg",
                        header = T, stringsAsFactors = F)
repro_cor <- repro_cor[ ,c("p1","p2","gcov_int")]

# Function to create matrix 

createCorrMat <- function (gcov_df) {
  # replace munged filenames with phenotype names
  gcov_df$p1 <- gsub(".*munged_sumstats/", "" , gcov_df$p1)
  gcov_df$p1 <- gsub("_for_ldsc.sumstats.gz", "" , gcov_df$p1) 
  
  gcov_df$p2 <- gsub(".*munged_sumstats/", "" , gcov_df$p2)
  gcov_df$p2 <- gsub("_for_ldsc.sumstats.gz", "" , gcov_df$p2) 
  
  phen <- unique(c(gcov_df$p1, gcov_df$p2))
  n <- length(phen)
  mat <- matrix(NA, n, n)                    # create matrix
  rownames(mat) <- colnames(mat) <- phen    # set col/rownames
  
  for (i in phen) {
    for (j in phen) {
      #################### FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!
      print(paste0(i, " x ", j))
      if (i == j) {
        int_val <- 1
      } else {
        int_val <- subset(gcov_df, p1==i & p2==j)$gcov_int
      }
      mat[i,j] <- int_val
    }
  }
  
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)] # sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
  mat <- round(cov2cor(mat), 5) # standardise
  return (mat)
}

# Apply to the hormone-infertility matrix
hi_mat <- createCorrMat(hi_cor)
write.table(hi_mat, 
            "/well/lindgren/samvida/hormones_infertility/lava_local/hormones_infertility_sample_overlap_from_ldsc.txt", 
            quote = F)   

# Apply to repro mat
repro_mat <- createCorrMat(repro_cor)
write.table(repro_mat, 
            "/well/lindgren/samvida/hormones_infertility/lava_local/repro_traits_sample_overlap_from_ldsc.txt", 
            quote = F)   

