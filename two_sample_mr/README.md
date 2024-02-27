Scripts to perform two-sample Mendelian randomisation (using summary statistics):

1. **1_grep_exposures_outcomes.sh** - Get summary statistics for the outcomes (infertility) at variants defined in the exposure instruments (obesity - BMI, WHR, and WHRadjBMI and hormones - FSH and testosterone). 
2. **2_harmonise.R** - Match effect alleles across summary statistics.
3. **3_perform_MR.R** - MR using the TwoSampleMR package in MRBase. Report results from different methods, test for heterogeneity, etc.
4. **4_plot_results.R** - Tables and plots in manuscript.
