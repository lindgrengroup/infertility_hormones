# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(tidyverse)
theme_set(theme_bw())

# Read in files ----

mainpath <- "/well/lindgren/samvida/hormones_infertility/exome_seq_results/conditioned_on_common"

HORMONES <- c("Oestradiol_F", "Testosterone_F",
              "Testosterone_M")

PTHRESH <- 1E-07

# Prepare columns for reading
to_numeric <- c("CHR", "POS", 
                "AC_Allele2", "AF_Allele2",
                "MissingRate", "BETA", "SE", "Tstat",
                "var", "p.value",
                "BETA_c", "SE_c", "Tstat_c", 
                "var_c", "p.value_c", "N")

# results
full_dat <- lapply(HORMONES, function (pheno) {
  print(pheno)
  # Read across all chromosomes
  fnames <- list.files(path = paste0(mainpath, "/results"),
                       pattern = paste0("chr.*_", pheno, "_.*"),
                       full.names = T)
  res <- lapply(fnames, function (fnm) {
    df <- read.table(fnm,
               sep = "\t", header = T, stringsAsFactors = F,
               colClasses = "character")
    # only retain the variants that were (unconditioned) exome-wide significant for this trait
    df <- df %>%
      mutate(across(all_of(to_numeric), as.numeric)) %>%
      filter(p.value <= PTHRESH)
    return (df)
  })
  res <- bind_rows(res)
  res <- res %>%
    mutate(phenotype = pheno) %>%
    arrange(CHR, POS)
  return (res)
})
full_dat <- bind_rows(full_dat)

# Write for manuscript table
to_write <- full_dat %>%
  mutate(uncond_beta = paste0(signif(BETA, 3), " (", signif(SE, 3), ")"),
         uncond_p = signif(p.value, 3),
         cond_beta = paste0(signif(BETA_c, 3), " (", signif(SE_c, 3), ")"),
         cond_p = signif(p.value_c, 3)) %>%
  select(all_of(c("phenotype", "MarkerID", 
                  "uncond_beta", "uncond_p",
                  "cond_beta", "cond_p")))

write.table(to_write, paste0(mainpath, "/results/for_manuscript_table.txt"),
            sep = "\t", quote = F, row.names = F)
