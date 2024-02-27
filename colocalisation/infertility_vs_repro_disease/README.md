Scripts must be run in the following order:

1. **1_prep_coloc_sumstats.sh** - For reproductive disease (endometriosis, PCOS, heavy menstrual bleeding, uterine fibroids) summary statistics on hg19, liftover to hg38 and add sample sizes. Create windows of +/- 50kb around each lead infertility variant to assess for colocalisation with a reproductive disease signal. Subset reproductive disease summary statistics to this window.
2. **2_run_coloc.R** - Munge sumstats from both sets of traits to match alleles and run colocalisation using the coloc package. Parallelise across every combination of infertility x reproductive disease.
3. **3_examine_coloc_results_repro_infertility.R** - Create table to report colocalisation results. 