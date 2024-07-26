Scripts in *prep_mvp_sumstats/*:

- **filter_plot_mvp_results.R** - Standardise the columns, ensure the X chromosome is coded the same in all studies, remove indels, remove SNPs with implausibly large standard errors (>10), remove duplicate entries, and print the summary statistics files formatted for meta-analysis in METAL. These scripts also output Manhattan plots and QQ plots (split by MAF bin). *Note that since the QQ plots show inflation at low MAF in the multi-ancestry meta-analysis, these files are also filtered to remove MAF less than 1%.
Submit with **submit_filter_plot_mvp_results.sh**
