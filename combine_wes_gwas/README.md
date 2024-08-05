Scripts to compare results from GWAS meta-analyses to WES analyses, to be run in the following order:

1. **1_extract_sig_hits.sh** - This needs to be run first to get rare variants, gene-based results, and common variants that pass their respective significance thresholds for each trait. 
2. **2_plot_combined_results.R** - Annotate the variant (both common and rare) results to their nearest genes; wrangle data to make the effect size vs MAF plots reported in the manuscript - by gene, by strata, etc.

In response to a reviewer request to replicate the rare variant analyses (from exomes in UK Biobank), the following scripts were added:

3. **3_check_wes_variants_in_gwas.sh** - Subset GWAS (all-ancestry and EUR) meta-analysis results for variants that reach significance in UKB exome sequencing. 
4. **4_make_manuscript_tables.R** - Flip the effect size to match minor allele and report summary statistics.