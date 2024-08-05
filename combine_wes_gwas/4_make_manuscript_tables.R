# Author: Samvida S. Venkatesh
# Date: 05/08/24

library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

HORMONES <- c("Oestradiol_F", "Testosterone_F", "Testosterone_M") 

rep_dat <- lapply(HORMONES, function (hr) {
  df <- read.table(paste0(hr, "_gwas_sumstats_wes_variants.txt"),
                   sep = "\t", header = T, stringsAsFactors = F)
  df$strata <- hr
  return (df)
})
names(rep_dat) <- HORMONES

# Process the data to:
# 1. align beta with minor allele
# 2. sort by chr:pos

proc_dat <- lapply(HORMONES, function (hr) {
  df <- rep_dat[[hr]]
  df <- df %>%
    mutate(flip_beta = MAF!=Freq1,
           NEW_BETA = ifelse(flip_beta, -BETA, BETA),
           write_beta = paste0(signif(NEW_BETA, 3), " (",
                               signif(SE, 3), ")"),
           write_maf = signif(MAF, 3),
           write_pval = signif(PVALUE, 3)) 
  
  df_for_manu <- df %>%
    arrange(CHROM, GENPOS) %>%
    select(all_of(c("ID", "write_maf", "write_beta", "write_pval",
                    "strata")))
  return (df_for_manu)
})

proc_dat <- bind_rows(proc_dat)
write.table(proc_dat, "for_manuscript.txt",
            sep = "\t", row.names = F, quote = F)
