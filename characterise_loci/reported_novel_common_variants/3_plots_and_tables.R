# Author: Samvida S. Venkatesh
# Date: 22/06/2022

library(tidyverse)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/conditional_analysis"

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone") 
SEX_STRATA <- c("F", "M", "sex_comb")
ANCESTRY_GROUPS <- c("EUR", "all_anc")

# Read classified lead SNPs ----

classified_lead_snps <- read.table(paste0(mainpath, 
                                          "/classified_all_lead_snps_all_strata.txt"),
                                   sep = "\t", header = T, stringsAsFactors = F)
# Make sure columns are in the right types
CHAR_COLS <- c("ID", "RSID", "Allele1", "Allele2", 
               "Direction", "VARID_1KG", "hormone", "sex_strata", "ancestry",
               "classification", "strata")
NUM_COLS <- c("CHROM", "GENPOS", "MAF", "Freq1", "FreqSE",
              "BETA", "SE", "PVALUE", "HetPVal")

classified_lead_snps <- classified_lead_snps %>%
  mutate(across(any_of(CHAR_COLS), as.character)) %>%
  mutate(across(any_of(NUM_COLS), as.numeric)) %>%
  mutate(sex_strata = recode(sex_strata, "FALSE" = "F"))

# Plot the number of reported and novel SNPs per strata ----

custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("novel_any_hormone", 
                                 "novel_this_hormone", "reported")
plotBars <- function (df) {
  summ_df <- df %>% 
    group_by(strata) %>%
    count(classification)
  
  resplot <- ggplot(summ_df, 
         aes(x = strata, y = n)) +
    geom_col(aes(fill = classification), width = 0.7) +
    scale_fill_manual(values = custom_three_diverge) +
    labs(x = "Strata", y = "Number of lead SNPs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  return (resplot)
}

lapply(HORMONES, function (hr) {
  for_plot <- classified_lead_snps %>% filter(hormone == hr & sex_strata != "sex_comb")
  
  # Barplots
  ggsave(paste0(mainpath, "/plots/", hr, "_number_snps_novel_reported_for_manuscript.png"),
         plotBars(for_plot), units = "in", height = 7, width = 7)
})

# Read annotated SNPs with functional predictions ----

annot_lead_snps <- read.table(paste0(mainpath, 
                                     "/annotated_classified_all_lead_snps_all_strata.txt"),
                              sep = "\t", header = T, stringsAsFactors = F)

# Retain only non-NA max fpred scores for each variant
# and create a note column with gene name (score)
annot_lead_snps <- annot_lead_snps %>%
  mutate(for_note = paste0(hgnc_symbol, " - ", fpred_max_label, " (", fpred_max_score, ")")) %>%
  group_by(ID, strata) %>%
  arrange(desc(fpred_max_score)) %>%
  mutate(note = str_c(for_note, collapse = "; ")) %>%
  select(-any_of(c("for_note", "hgnc_symbol", "gene_id",
                   "fpred_max_score", "fpred_max_label"))) %>%
  mutate(note = gsub(" - NA \\(NA\\)", "", note)) %>%
  distinct()

# Bind in rows for the SNPs which are not annotated
not_annot <- classified_lead_snps %>%
  filter(!ID %in% annot_lead_snps$ID)
annot_lead_snps <- bind_rows(annot_lead_snps, not_annot)

# Annotate SNPs with other associated traits ----

pub_variants <- read.table(paste0(mainpath, "/published_snps/all_hormones_reported.txt"),
                                             sep = " ", header = F, stringsAsFactors = F)
colnames(pub_variants) <- c("ID", "VARID_1KG", "reported_trait")

pub_variants <- pub_variants %>% 
  group_by(ID) %>% 
  summarize(reported_trait = paste0(reported_trait, collapse = "; "))

annot_lead_snps <- annot_lead_snps %>%
  left_join(pub_variants, by = "ID") %>%
  group_by(strata) %>%
  arrange(CHROM, GENPOS, by_group = T)

write.table(annot_lead_snps,
            paste0(mainpath, "/with_gwascat_opentargets_annotated_classified_lead_snps_all_strata.txt"), 
            sep = "\t", row.names = F, quote = F)


