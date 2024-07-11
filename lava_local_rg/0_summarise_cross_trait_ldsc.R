# Author: Samvida S. Venkatesh
# Date: 08/07/2024
# Adapted from: https://github.com/josefin-werme/LAVA/blob/main/vignettes/sample_overlap.Rmd

library(tidyverse)

mainpath <- "/well/lindgren/samvida/hormones_infertility/lava_local"

# Read data ----

# Read mapping table for filename to pheno_name

filename_map <- read.table(paste0(mainpath, "/input_sumstats.txt"),
                           sep = "\t", header = T, stringsAsFactors = F)
# clean phenotype names
# replace special characters with "_"
filename_map$phenotype <- gsub("[^A-Za-z0-9_]", "_", filename_map$phenotype)
# replace any double-underscores with "_"
filename_map$phenotype <- gsub("__", "_", filename_map$phenotype)
# replace uppercase with lowercase
filename_map$phenotype <- tolower(filename_map$phenotype)

NTARGETS <- nrow(filename_map) - 1

# Read in tables of LDSC cross-trait rG results

ldsc_results <- lapply(0:NTARGETS, function (target) {
  res <- read.table(paste0(mainpath, "/cross_trait_ldsc/", target, "_vs_all_phenos_rg.r2"),
                    header = T, stringsAsFactors = F)
  res <- res[, c("p1", "p2", "gcov_int")]
  
  # Replace filenames with phenotype names
  res$p1 <- filename_map$phenotype[match(res$p1, filename_map$filename)]
  res$p2 <- filename_map$phenotype[match(res$p2, filename_map$filename)]
  return (res)
})
ldsc_results <- bind_rows(ldsc_results)

# Function to create matrix ----

createCorrMat <- function (gcov_df) {
  phen <- unique(c(gcov_df$p1, gcov_df$p2))
  phen <- phen[phen %in% PHENOS]
  n <- length(phen)
  mat <- matrix(NA, n, n)                    # create matrix
  rownames(mat) <- colnames(mat) <- phen    # set col/rownames
  
  for (i in phen) {
    for (j in phen) {
      print(paste0(i, " x ", j))
      if (i == j) {
        int_val <- 1
      } else {
        int_val <- subset(gcov_df, p1==i & p2==j)$gcov_int
        if (length(int_val) == 0) {
          int_val <- subset(gcov_df, p1==j & p2==i)$gcov_int
        }
      }
      mat[i,j] <- int_val
    }
  }
  mat <- round(cov2cor(mat), 5) # standardise
  return (mat)
}

# Apply to the cross-trait results
ldsc_mat <- createCorrMat(ldsc_results)
write.table(ldsc_mat, 
            "/well/lindgren/samvida/hormones_infertility/lava_local/cross_trait_sample_overlap_from_ldsc.txt", 
            quote = F)   

