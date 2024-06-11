library(tidyverse)

# Read coloc results ----

REPRO_STRATA <- c("endometriosis", "PCOS")
INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis2",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5")

coloc_res <- read.table("infertility_repro_coloc.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

coloc_res <- coloc_res %>%
  mutate(ml_ratio = PP.H4/PP.H3) %>%
  select(all_of(c("lead_variant", "nsnps", "PP.H3", "PP.H4", "ml_ratio",
                  "infertility", "repro_trait")))

to_write <- coloc_res %>%
  mutate(PP.H3 = paste0(signif(PP.H3*100, 3), "%"),
         PP.H4 = paste0(signif(PP.H4*100, 3), "%"),
         ml_ratio = signif(ml_ratio, 3)) %>%
  filter(repro_trait %in% REPRO_STRATA) %>%
  pivot_wider(id_cols = c("lead_variant", "infertility"),
              names_from = repro_trait,
              values_from = c("nsnps", "PP.H3", "PP.H4", "ml_ratio")) %>%
  arrange(infertility, lead_variant) %>%
  select(all_of(c("lead_variant", "infertility",
                  "nsnps_endometriosis", "PP.H3_endometriosis", "PP.H4_endometriosis", "ml_ratio_endometriosis",
                  "nsnps_PCOS", "PP.H3_PCOS", "PP.H4_PCOS", "ml_ratio_PCOS")))
  
write.table(to_write, "for_thesis_table_repro_infertility_coloc_results.txt",
            sep = "\t", row.names = F, quote = F)
