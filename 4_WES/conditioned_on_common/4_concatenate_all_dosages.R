# Author: Samvida S. Venkatesh
# Date: 20/06/24

library(tidyverse)

mainpath <- "/well/lindgren/samvida/hormones_infertility/exome_seq_results/conditioned_on_common"

lead_variants <- read.table(paste0(mainpath, "/nearby_common_variants/unique_rsids_to_get.txt"),
                            sep = "\t", header = F, stringsAsFactors = F)
colnames(lead_variants) <- c("SNP", "CHR", "POS")

dosage_files <- lapply(lead_variants$SNP, function (rsid) {
  df <- read.table(paste0(mainpath, "/dosage_files/", rsid, "_dosages.txt"),
                   sep = " ", header = T, stringsAsFactors = F)
  df <- df[, c(1,2,5)]
  
  # Threshold genotypes
  dosages <- df[,3]
  geno <- ifelse(dosages > 0 & dosages < 0.5, 0,
                 ifelse(dosages > 0.5 & dosages < 1, 1, 
                        ifelse(dosages > 1 & dosages < 1.5, 1,
                           ifelse(dosages > 1.5 & dosages < 2, 2, dosages))))
  
  res <- data.frame(FID = df[,1],
                    IID = df[, 2],
                    genotype = geno)
  colnames(res) <- c("FID", "IID", rsid)
  return (res)
})

dosage_covars <- reduce(dosage_files, 
                        function(x, y) inner_join(x, y, by = c("FID", "IID")))

write.table(dosage_covars, 
            paste0(mainpath, "/snp_dosages_to_condition.txt"), sep = "\t",
            quote = F, row.names = F)
