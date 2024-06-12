Added in response to reviewer request, look at the INHBB locus for related traits: AMH, PCOS, age at menopause, and F-ANOV

Scripts must be run in the following order:

1. **1_prep_coloc_inhbb_sumstats.sh** - For Anti-Mullerian hormone summary statistics on hg19 (downloaded from Verdiesen et al. 2023), liftover to hg38 and add sample sizes. Create windows of +/- 50kb around the INHBB lead infertility variant to assess for colocalisation with an AMH signal. Subset AMH summary statistics to this window.
2. **2_run_coloc_inhbb.R** - Munge sumstats from all traits (PCOS, AMH, age at menopause, and F-ANOV) to match alleles and run colocalisation using the coloc package pairwise across all combinations. 
3. **3_examine_inhbb.R** - Create figure to report results.