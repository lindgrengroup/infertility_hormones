Scripts in this folder:

Scripts need to be run in the following order:

1. **1_get_data.R** - Extract codes for all hormones of interest from full primary care data linked to UK Biobank, and perform initial QC to retain only biomarkers with > 1 measurement in > 1000 individuals in primary care (for longitudinal resource construction). 
2. **2_popn_QC.R** - Perform population-level QC to exclude records based on the following filters:
	- age not between 20-80 yrs
	- implausible values defined as outside +/-10% of the UKBB assessment centre minimum and maximum values
	- extreme values defined as > 5 SD away from population mean
3. **3_add_UKB.R** - Integrate UKBB assessment centre or biomarker measurements with linked primary care data (QC'd above) and filter to only retain individuals with longitudinal records (> 1 measurement of each biomarker).
4. **4_prep_cross_sec_dat.R** - Remove zero or negative values and any records that are extreme on the individual level (i.e. records which cause large jumps in the time-series) by taking the log-fold-change between consecutive measurements, adjusted for time between consecutive measurements; any observations that are involved in an extreme jump > 3 SD away from population mean jump are removed. Retain only the median value per individual (if there are multiple measurements). To prepare for GWAS, within men and women separately as well as sex-combined, adjust for covariates age, age-square, and data provider, and rank-based inverse normally transform the data. 

Scripts in *prep/*:

- **gp_age_sex_BMI_annotation.R** - Link primary care data (GP data) with UK Biobank main phenotypes file (sex and average BMI across multiple assessment centre visits). Also integrate approximate date of birth (from month and year of birth in main file) to calculate age at event from the date at event field in GP data.

Code-lists for hormones are in the folder *code_lists/*, including:

- Lists of biomarkers in **traits_available.txt** and **qcd_traits_available.txt**
- Corresponding primary-care read codes in *primary_care/* (adapted from [Denaxas et al.](https://github.com/spiros/ukb-biomarker-phenotypes)),
- Corresponding UKBB assessment centre codes in *main_phenotypes/* (manually curated)
- List of codes for bariatric surgery records (and history of) in **bariatric_surgery_codes.txt**

Files required for quality control are in the folder *qc/*, including:
- Minimum and maximum values of biomarkers available from UKB main assessment centre measurements in **ukb_min_max.txt**

In response to reviewers requests, the following scripts were added:
- **check_age_menopause.R** - to check the age at menopause distributions in women with/without FSH & LH measurements so we know how worried we should be about ascertainment.

- **4_prep_cross_sec_dat_normal_menopause.R** - Same data extraction steps, for FSH and LH only, for those women who report age at menopause >45 years so we can RINT this data separately. 