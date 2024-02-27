Scripts in this folder are to assess results from whole exome sequencing (WES) analyses in the UK Biobank:

1. **1_prep_traits_for_wes.R** - For quantitative traits, creates a single file with values for all traits and covariates. 
2. **2_filter_plot_[gene/variant]_results.R** - QQ and Manhattan plots for gene-based and variant-level tests within each MAF and annotation group. Submit with **2_submit_filter_plot_results.sh** and batch submit across all strata with **batch_submit_filter_plot_results.sh**

- **concatenate_results.sh** - merge the per-chromosome files into a single file per trait. 
- **subset_gene_results_significant.R** - to create a table showing the different test (MAF threshold and annotation category) results per gene, for any gene that has a significant (P<5E-06) association with the trait.

Scripts to perform variant QC, sample QC, determine superpopulation labels, relatedness, and association analyses for WES tests are all documented in Duncan Palmer's repository here: https://github.com/astheeggeggs/BRaVa_curation.
