library(tidyverse)

# Read coloc results ----

INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5",
                   "male_infertility")

map_hypothesis_desc <- data.frame(hypothesis = c("H0", "H1", "H2", "H3", "H4"),
                                  hypo_desc = c("no_assoc", 
                                                "only_infert", "only_gene", 
                                                "both_diff_causal", "both_same_causal"))

coloc_res <- lapply(INFERT_STRATA, function (infert) {
  res <- read.table(paste0(infert, "_eur_gtex_coloc_results.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  return (res)
})
coloc_res <- bind_rows(coloc_res)

coloc_res <- coloc_res %>%
  mutate(ml_ratio = PP.H4.abf/PP.H3.abf)

coloc_res_sig <- coloc_res %>%
  filter(PP.H4.abf >= 0.5) 

to_write <- coloc_res_sig %>%
  mutate(chrpos_write = paste0("chr", Chromosome.scaffold.name, ":", 
                               Gene.start..bp., "-", 
                               Gene.end..bp.),
         PP_H3 = paste0(signif(PP.H3.abf*100, 3), "%"),
         PP_H4 = paste0(signif(PP.H4.abf*100, 3), "%"),
         ml_ratio = signif(ml_ratio, 3)) %>%
  select(all_of(c("chrpos_write", "Gene.name", "Tissue", "Trait",
                  "PP_H3", "PP_H4", "ml_ratio")))

write.table(to_write, "for_thesis_table_gtex_coloc_results.txt",
            sep = "\t", row.names = F, quote = F)
