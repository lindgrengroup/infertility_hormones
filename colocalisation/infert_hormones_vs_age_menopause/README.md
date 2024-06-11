Added in response to reviewer request

Scripts must be run in the following order:

1. **1_prep_coloc_sumstats.sh** - For age at menopause summary statistics on hg19 (downloaded from Ojavee et al. 2023), liftover to hg38 and add sample sizes. Create windows of +/- 50kb around each lead infertility variant to assess for colocalisation with an age-at-menopause signal. Subset age at menopause summary statistics to this window.
2. **2_run_coloc.R** - Munge sumstats from both sets of traits to match alleles and run colocalisation using the coloc package. Parallelise across every combination of infertility x age at menopause.
3. **3_examine_coloc_results_age_at_menopause.R** - Create table to report colocalisation results. 