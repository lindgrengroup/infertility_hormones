library(tidyverse)

# Read coloc results ----

HORMONE_STRATA <- c("FSH_F", "LH_F", "Progesterone_F",
                    "Oestradiol_F", "Testosterone_F")
INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis2",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5")

map_hypothesis_desc <- data.frame(hypothesis = c("H0", "H1", "H2", "H3", "H4"),
                                  hypo_desc = c("no_assoc", 
                                                "only_infert", "only_hormone", 
                                                "both_diff_causal", "both_same_causal"))

hr_infert_res <- lapply(INFERT_STRATA, function (infert) {
  per_hr <- lapply(HORMONE_STRATA, function (hr) {
    coloc_res <- read.table(paste0("coloc_results_", hr, "_EUR_", infert, "_EUR.txt"),
                            sep = "\t", header = T, stringsAsFactors = F)
    return (coloc_res)
  })
  per_hr <- bind_rows(per_hr)
  return (per_hr)
})
hr_infert_res <- bind_rows(hr_infert_res)

hr_infert_res <- hr_infert_res %>%
  mutate(ml_ratio = PP.H4.abf/PP.H3.abf)

hr_infert_res_sig <- hr_infert_res %>%
  filter(PP.H4.abf >= 0.5) 

to_write <- hr_infert_res_sig %>%
  mutate(chrpos_write = paste0("chr", chr_sentinel, ":", pos_sentinel - 50000,
                               "-", pos_sentinel + 50000),
         PP_H3 = paste0(signif(PP.H3.abf*100, 3), "%"),
         PP_H4 = paste0(signif(PP.H4.abf*100, 3), "%"),
         ml_ratio = signif(ml_ratio, 3)) %>%
  select(all_of(c("chrpos_write", "hormone", "infertility", "nsnps",
                  "PP_H3", "PP_H4", "ml_ratio")))

write.table(to_write, "for_thesis_table_coloc_results.txt",
            sep = "\t", row.names = F, quote = F)



