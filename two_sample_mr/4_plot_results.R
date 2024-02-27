# Author: Samvida S. Venkatesh
# Date: 27/06/23

library(tidyverse)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/two_sample_mr"

FEMALE_OUTCOMES <- c("female_infertility_analysis1", "female_infertility_analysis2", 
                     "female_infertility_analysis3", "female_infertility_analysis4", 
                     "female_infertility_analysis5")
MALE_OUTCOMES <- "male_infertility"

mr_res <- read.table(paste0(mainpath, "/results_230627.txt"),
            sep = "\t", header = T, stringsAsFactors = F)

# Wrangle results ----

# Add OR and 95% C.I.
mr_res_or <- mr_res %>%
  mutate(OR = exp(b), lci = exp(b - 1.96*se), uci = exp(b + 1.96*se),
         fdr_sig = ifelse(pval.fdr <= 0.05, "fdr_sig", "ns"),
         method = recode(method, "MR Egger" = "Egger",
                         "Weighted median" = "Wt med",
                         "Inverse variance weighted" = "IVW",
                         "Weighted mode" = "Wt mode"))

# Plot results ----

colpal_female <- c("#D35C79", "#666666")
names(colpal_female) <- c("fdr_sig", "ns")

colpal_male <- c("#009593", "#666666")
names(colpal_male) <- c("fdr_sig", "ns")

getORPlot <- function (mr_df, colpal_use = colpal_female) {
  or_plot <- ggplot(mr_df, 
                    aes(x = OR, y = method,
                        colour = fdr_sig, fill = fdr_sig)) +
    facet_wrap(~exposure, ncol = 2, scales = "free_x") +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        colour = fdr_sig, fill = fdr_sig),
                    size = 0.7,
                    position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = colpal_use, guide = "none") +
    scale_colour_manual(values = colpal_use, guide = "none") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) 
  return (or_plot)
}

lapply(FEMALE_OUTCOMES, function (fout) {
  df_plot <- mr_res_or %>% filter(outcome == fout)
  ggsave(paste0("plots/outcome_", fout, ".png"), 
         getORPlot(df_plot, colpal_use = colpal_female), 
         height = 7, width = 7, units = "in")
})

male_plot <- mr_res_or %>% filter(outcome == "male_infertility")
ggsave("plots/outcome_male_infertility.png", 
       getORPlot(male_plot, colpal_use = colpal_male), 
       height = 7, width = 7, units = "in")

