# Author: Samvida S. Venkatesh
# Date: 19/06/24

library(tidyverse)

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone")
SEX_STRATA <- c("F", "M")
ANC_GROUPS <- c("EUR", "all_anc")

mainpath <- "/well/lindgren/samvida/hormones_infertility/meta_results_230613/lead_snps"

NUMERIC_COLS <- c("CHROM", "GENPOS", "MAF", "Freq1", "FreqSE", 
                  "BETA", "SE", "PVALUE", "HetPVal")
CHAR_COLS <- c("ID", "RSID", "Allele1", "Allele2", "Direction",
               "hormone", "sex_strata", "ancestry")

all_dat <- lapply(HORMONES, function (hr) {
  res_list <- lapply(SEX_STRATA, function (ss) {
    res_df <- lapply(ANC_GROUPS, function (anc) {
      df <- NULL
      fname <- paste0(mainpath, "/all_lead_snp_sumstats_",
                      hr, "_", ss, "_", anc, "_with_rsids.txt")
      if (file.exists(fname)) {
        df <- read.table(fname,
                         sep = "\t", header = T, stringsAsFactors = F) %>%
          mutate(across(all_of(NUMERIC_COLS), as.numeric),
                 across(all_of(CHAR_COLS), as.character))
      } 
      return (df)
    })
    res_df <- bind_rows(res_df)
    return (res_df)
  })
  res_list <- bind_rows(res_list)
  return (res_list)
})
all_dat <- bind_rows(all_dat)

perc_het <- length(which(all_dat$HetPVal < 0.05))/nrow(all_dat)

print(paste0("Percent SNPs not heterogeneous: ", signif((1-perc_het)*100,3), "%"))
