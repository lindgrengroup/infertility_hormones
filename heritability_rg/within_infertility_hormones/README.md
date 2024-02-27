Scripts to be run in the following order:

1. **1_munge_sumstats_ldsc.sh** - Prepare summary statistics for heritability and rG analyses using LDSC. 
2. **2_calculate_rg_h2.sh** - Use LDSC to calculate rG between all infertility and hormone traits, as this also reports the heritability. Only report rG for pairs of traits with significant heritability (Z>4). 
3. **3_plot_heritability_rg.R** - Create the plots in the manuscript.

These additional scripts were used for calculations reported in the paper:

- **check_ukb_all_pheno_heritability.R** - To estimate how low the heritability of infertility is compared to other diseases of similar prevalence in UK Biobank, we looked at the Neale lab estimates of h2 for all binary traits with prevalence <5%. 
- **convert_observed_liability_scale_h2.R** - Convert the binary (infertility) heritability estimate from observed scale produced by LDSC to liability scale assuming prevalence in the meta-analysis here is the same as population prevalence (as there are no good estimates of population prevalence, especially for various subtypes of female infertility).
