Scripts in this folder to perform GWASs in UK Biobank for the quantitative hormone traits:

1. **1_sample_QC.R** and **1_sample_QC_non_eur.R**- Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, samples with poor heterozygosity or missingness, etc.). Save sample ids that pass QC along with genotyping-related covariates, i.e. genotyping array and UKB assessment centre.
2. **2_submit_regenie_jobs.R** - Loops over all hormones and sex-specific analyses strata to submit REGENIE GWAS. Relies on helper scripts in *./regenie_scripts/*. 

To address reviewer comments, we added the following scripts to subset women with age at menopause >45 years for separate GWAS.
1. **1a_sample_QC_normal_menopause.R** - Re-run sample QC for women with normal age at menopause (age>45 yrs)

Scripts in *./regenie_scripts/*:

1. **regenie_step1.sh** - For REGENIE null model fitting in step 1.
2. **regenie_step2.sh** and **regenie_step2_xchr.sh** - For REGENIE model fitting in step 2, array job split by chromosome. X chromosome requires a different script because the UKB sample file is not the same as the single sample file for the autosomes.


