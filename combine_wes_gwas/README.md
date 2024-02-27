Scripts to compare results from GWAS meta-analyses to WES analyses, to be run in the following order:

1. **1_extract_sig_hits.sh** - This needs to be run first to get rare variants, gene-based results, and common variants that pass their respective significance thresholds for each trait. 
2. **2_plot_combined_results.R** - Annotate the variant (both common and rare) results to their nearest genes; wrangle data to make the effect size vs MAF plots reported in the manuscript - by gene, by strata, etc.

