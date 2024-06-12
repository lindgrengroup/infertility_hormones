library(tidyverse)

# Read coloc results ----

INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis2",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5")

coloc_res <- read.table("age_at_menopause_infertility_coloc_results.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

coloc_res <- coloc_res %>%
  mutate(ml_ratio = PP.H4/PP.H3) %>%
  select(all_of(c("lead_variant", "nsnps", "PP.H3", "PP.H4", "ml_ratio",
                  "infertility")))

to_write <- coloc_res %>%
  mutate(PP.H3 = paste0(signif(PP.H3*100, 3), "%"),
         PP.H4 = paste0(signif(PP.H4*100, 3), "%"),
         ml_ratio = signif(ml_ratio, 3)) %>%
  arrange(infertility, lead_variant) %>%
  select(all_of(c("lead_variant", "infertility",
                  "nsnps", "PP.H3", "PP.H4", "ml_ratio")))
  
write.table(to_write, "manuscript_table_aam_infertility_coloc_results.txt",
            sep = "\t", row.names = F, quote = F)
