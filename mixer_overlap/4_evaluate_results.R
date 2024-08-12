# Author: Samvida S. Venkatesh
# Date: 09/0824

library(tidyverse)

# Read MiXeR results ----

PHENOS <- read.table("input_sumstats.txt", sep = "\t",
                     header = T, stringsAsFactors = F)$phenotype

# Univariate results from .test ----

univar_res <- lapply(c(1:3), function (i) {
  per_target <- lapply(c(4:length(PHENOS)), function (j) {
    df <- read.table(paste0("results/", PHENOS[i], "_and_", PHENOS[j], ".test.csv"),
                     sep = "\t", header = T, stringsAsFactors = F)
    df <- df %>%
      mutate(pheno = gsub("results/", "", fname),
             pheno = gsub("\\.test\\.json", "", pheno)) 
    # replace the AIC/BIC values with those from .fit
    df_fit <- read.table(paste0("results/", PHENOS[i], "_and_", PHENOS[j], ".fit.csv"),
                                     sep = "\t", header = T, stringsAsFactors = F)
    df_fit <- df_fit %>%
      mutate(pheno = gsub("results/", "", fname),
             pheno = gsub("\\.fit\\.json", "", pheno)) 

    df$AIC <- df_fit$AIC[match(df$pheno, df_fit$pheno)]
    df$BIC <- df_fit$BIC[match(df$pheno, df_fit$pheno)]
    return (df)
  })
  per_target <- bind_rows(per_target)
  return (per_target)
})
univar_res <- bind_rows(univar_res)

univar_res_summ <- univar_res %>%
  group_by(pheno) %>%
  summarise(across(-any_of(c("fname", "pheno")), mean))

write_univar_res <- univar_res_summ %>%
  mutate(pi = paste0(signif(pi..mean.*100, 3), "% (", signif(pi..std.*100, 3), ")"),
         sig2_beta = paste0(signif(sig2_beta..mean., 3), " (", signif(sig2_beta..std., 3), ")"),
         sig2_zero = paste0(signif(sig2_zero..mean., 3), " (", signif(sig2_zero..std., 3), ")"),
         n_causal = paste0(signif(nc.p9..mean., 3), " (", signif(nc.p9..std., 3), ")"),
         h2 = paste0(signif(h2..mean.*100, 3), "% (", signif(h2..std.*100, 3), ")"),
         AIC = signif(AIC, 3), BIC = signif(BIC, 3),
         pheno = factor(pheno, levels = PHENOS)) %>%
  arrange(pheno) %>%
  select(all_of(c("pheno", "pi", "sig2_beta", "sig2_zero", 
                  "n_causal", "h2", "AIC", "BIC")))

write.table(write_univar_res, "for_table_univariate_results.txt",
            sep = "\t", row.names = F, quote = F)

# Bivariate results ----

bivar_res <- lapply(c(1:3), function (i) {
  per_target <- lapply(c(4:length(PHENOS)), function (j) {
    df <- read.table(paste0("results/", PHENOS[i], "_vs_", PHENOS[j], ".csv"),
                     sep = "\t", header = T, stringsAsFactors = F)
    # keep AIC/BIC from fit, rest from test
    df[2, c("best_vs_min_AIC", "best_vs_min_BIC", 
            "best_vs_max_AIC", "best_vs_max_BIC")] <- df[1, c("best_vs_min_AIC", "best_vs_min_BIC", 
                                                              "best_vs_max_AIC", "best_vs_max_BIC")]
    return (df[2, ])
  })
  per_target <- bind_rows(per_target)
  return (per_target)
})
bivar_res <- bind_rows(bivar_res) %>%
  mutate(pheno = gsub("results/", "", fname),
         pheno = gsub("\\.test\\.json", "", pheno),
         trait1 = gsub("_vs_.*", "", pheno),
         trait2 = gsub(".*_vs_", "", pheno)) 

write_bivar_res <- bivar_res %>%
  mutate(nc12 = paste0(signif(nc12.p9..mean., 3), " (", signif(nc12.p9..std., 3), ")"),
         nc1 = paste0(signif(nc1.p9..mean., 3), " (", signif(nc1.p9..std., 3), ")"),
         nc2 = paste0(signif(nc2.p9..mean., 3), " (", signif(nc2.p9..std., 3), ")"),
         rho_beta = paste0(signif(rho_beta..mean., 3), " (", signif(rho_beta..std., 3), ")"),
         rg = paste0(signif(rg..mean., 3), " (", signif(rg..std., 3), ")"),
         best_vs_min_AIC = signif(best_vs_min_AIC, 3), 
         best_vs_min_BIC = signif(best_vs_min_BIC, 3),
         best_vs_max_AIC = signif(best_vs_max_AIC, 3),
         best_vs_max_BIC = signif(best_vs_max_BIC, 3),
         trait1 = factor(trait1, levels = PHENOS),
         trait2 = factor(trait2, levels = PHENOS)) %>%
  arrange(trait1, trait2) %>%
  select(all_of(c("trait1", "trait2",
                  "nc12", "nc1", "nc2", 
                  "rho_beta", "rg",
                  "best_vs_min_AIC", "best_vs_min_BIC", "best_vs_max_AIC", "best_vs_max_BIC")))

write.table(write_bivar_res, "for_table_bivariate_results.txt",
            sep = "\t", row.names = F, quote = F)

# For interpretation ----

comp_ldsc <- read.table("bivariate_results_with_ldsc_rg.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

comp_ldsc$rg_mixer <- as.numeric(gsub("\\(.*", "", comp_ldsc$rg))
comp_ldsc$se_mixer <- as.numeric(gsub(")", "", gsub(".*\\(.", "", comp_ldsc$rg)))
comp_ldsc$rg_ldsc <- as.numeric(gsub("\\(.*", "", comp_ldsc$ldsc_rg))
comp_ldsc$se_ldsc <- as.numeric(gsub(")", "", gsub(".*\\(.", "", comp_ldsc$ldsc_rg)))

comp_ldsc <- comp_ldsc %>%
  mutate(het_z = (rg_mixer-rg_ldsc)/sqrt(se_mixer^2 + se_ldsc^2),
         rg_p = 2*pnorm(abs(het_z), lower.tail = F),
         p_write = signif(rg_p, 3))

write.table(comp_ldsc, "bivariate_results_with_ldsc_rg.txt",
            sep = "\t", row.names = F, quote = F)

bivar_annot <- bivar_res %>%
  mutate(prop_shared = nc12.p9..mean. / (nc12.p9..mean. + nc1.p9..mean. + nc2.p9..mean.),
         prop_unique_other = nc2.p9..mean. / (nc12.p9..mean. + nc1.p9..mean. + nc2.p9..mean.),
         prop_unique_infert = nc1.p9..mean. / (nc12.p9..mean. + nc1.p9..mean. + nc2.p9..mean.),
         rg_p = 2*pnorm(abs(rg..mean./rg..std.), lower.tail = F),
         shared_rg_p = 2*pnorm(abs(rho_beta..mean./rho_beta..std.), lower.tail = F)) 

##### TO DO EVALUATE THESE FOR SUPP TEXT

# high polygenic overlap: prop_shared > 25%
high_overlap <- bivar_annot %>% filter(prop_shared > 0.25)

# low polygenic overlap (<25%) but high local rG (P<0.05/42)
locally_specific <- bivar_annot %>% filter(prop_shared <= 0.25 & shared_rg_p < 0.05) 

 
