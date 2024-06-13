# Author: Samvida S. Venkatesh
# Date: 17/10/23

library(tidyverse)
library(coloc)

TRAITS <- c("female_infertility_analysis3",
            "PCOS",
            "age_at_menopause",
            "AMH")

mainpath <- "/well/lindgren/samvida/hormones_infertility/colocalisation"

# Read data ----

inhbb_snp <- "chr2_120388925"

# Load GWAS results 

infert <- read.table(paste0(mainpath, "/sentinel_windows/", inhbb_snp,
                            "_female_infertility_analysis3.txt"),
                     sep = "\t", header = T, stringsAsFactors = F)
infert$Allele1 <- toupper(infert$Allele1)
infert$Allele2 <- toupper(infert$Allele2)
infert$markerMatch <- paste0("chr", infert$chr, "_", infert$pos)

pcos <- read.table(paste0(mainpath, "/sentinel_windows/", inhbb_snp,
                          "_PCOS.txt"),
                   sep = "\t", header = T, stringsAsFactors = F)
pcos$markerMatch <- paste0(pcos$CHR, "_", pcos$POS)

aam <- read.table(paste0(mainpath, "/sentinel_windows/", inhbb_snp, "_age_at_menopause.txt"), 
                  sep = "\t", header = T, stringsAsFactors = F)
aam$markerMatch <- paste0(aam$CHR, "_", aam$POS)

amh <- read.table(paste0(mainpath, "/sentinel_windows/", inhbb_snp, "_verdiesen_AMH.txt"), 
                  sep = "\t", header = T, stringsAsFactors = F)
amh$markerMatch <- paste0(amh$CHR, "_", amh$POS)

# Get AMH vs infertility colocalisation results ----

# Wrangle infertility results retaining only common >1% markers
inf_df <- infert %>% filter(maf >= 0.01) %>%
  mutate(varbeta_infert = StdErr^2,
         N_infert = N_CASES + N_CONTROLS,
         casefrac_infert = N_CASES/N_infert) %>%
  rename(MarkerName_infert = MarkerName,
         A1_infert = Allele1, A2_infert = Allele2,
         MAF_infert = maf, BETA_infert = Effect, 
         PVAL_infert = P.value) %>%
  select(all_of(c("MarkerName_infert", "A1_infert", "A2_infert",
                  "MAF_infert", "BETA_infert", "varbeta_infert",
                  "PVAL_infert", "casefrac_infert", "N_infert",
                  "markerMatch")))

# Wrangle AAM results retaining only common >1% markers
amh_df <- amh %>% 
  mutate(MAF = ifelse(Effect_Allele_Frequency > 0.5, 1-Effect_Allele_Frequency, Effect_Allele_Frequency),
         varbeta_amh = SE^2,
         N_amh = N) %>%
  filter(MAF >= 0.01) %>%
  rename(MarkerName_amh = RSID,
         A1_amh = Effect_Allele, A2_amh = Other_Allele,
         MAF_amh = MAF, BETA_amh = BETA, 
         PVAL_amh = PVALUE) %>%
  select(all_of(c("MarkerName_amh", "A1_amh", "A2_amh",
                  "MAF_amh", "BETA_amh", "varbeta_amh", 
                  "PVAL_amh", "N_amh", 
                  "markerMatch")))

# Merge 
full_df <- inner_join(inf_df, amh_df, by = "markerMatch")

d_infert <- list(snp = full_df$markerMatch, 
                 pvalues = full_df$PVAL_infert, 
                 beta = full_df$BETA_infert, 
                 varbeta = full_df$varbeta_infert, 
                 MAF = full_df$MAF_infert, 
                 s = full_df$casefrac_infert,
                 N = full_df$N_infert, type = "cc")

d_amh <- list(snp = full_df$markerMatch, 
              pvalues = full_df$PVAL_amh, 
              beta = full_df$BETA_amh, 
              varbeta = full_df$varbeta_amh, 
              MAF = full_df$MAF_amh,
              N = full_df$N_amh, type = "quant")

coloc_df <- coloc.abf(dataset1 = d_infert, 
                      dataset2 = d_amh, p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)

# Extract results
amh_infert_coloc_res <- data.frame(nsnps = coloc_df$summary[1],
                                   PP.H0 = coloc_df$summary[2],
                                   PP.H1 = coloc_df$summary[3],
                                   PP.H2 = coloc_df$summary[4],
                                   PP.H3 = coloc_df$summary[5],
                                   PP.H4 = coloc_df$summary[6])

# Forest plot for F-ANOV SNP with PCOS and age at menopause ----

# effect allele = T
infert_plot <- infert %>% filter(markerMatch == "chr2_120388925") %>%
  rename(all_of(c(BETA = "Effect", SE = "StdErr", PVALUE = "P.value"))) %>%
  select(all_of(c("BETA", "SE", "PVALUE"))) %>%
  mutate(trait = "infertility")

# effect allele = T
pcos_plot <- pcos %>% filter(markerMatch == "chr2_120388925") %>%
  select(all_of(c("BETA", "SE", "PVALUE"))) %>%
  mutate(trait = "PCOS")

# effect allele = C
aam_plot <- aam %>% filter(markerMatch == "chr2_120388925") %>%
  select(all_of(c("BETA", "SE", "PVALUE"))) %>%
  mutate(BETA = -BETA, 
         trait = "age_at_menopause")

df_plot <- bind_rows(infert_plot, pcos_plot, aam_plot) %>%
  mutate(signif = ifelse(PVALUE <= 5E-08, "gws_sig", "ns"),
         lci = BETA - 1.96*SE, uci = BETA + 1.96*SE)

# Pretty plot
theme_set(theme_bw())
colpal_female <- c("#D35C79", "#666666")
names(colpal_female) <- c("gws_sig", "ns")

forest_plot <- ggplot(df_plot, 
                      aes(x = BETA, y = trait,
                          colour = signif, fill = signif)) +
  geom_pointrange(aes(xmin = lci, xmax = uci,
                      colour = signif, fill = signif),
                  size = 0.7,
                  position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = colpal_female, guide = "none") +
  scale_colour_manual(values = colpal_female, guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) 

ggsave(paste0(mainpath, "/inhbb_forest_plot.png"), 
       forest_plot, 
       height = 3, width = 3, units = "in")

