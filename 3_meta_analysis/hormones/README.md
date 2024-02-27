Before meta-analysis, summary statistics must be QCd.

Scripts in *filter_individual_results/*:

The scripts titled **filter_plot_indiv_results_[BOLT/deCode/regenie].R** are to standardise the columns from individual cohort summary statistics, ensure the X chromosome is coded the same in all studies, remove indels, remove SNPs with implausibly large standard errors (>10), remove duplicate entries, and print the summary statistics files formatted for meta-analysis in METAL. These scripts also output Manhattan plots and QQ plots (split by MAF bin).

The above scripts can be submitted individually with **submit_filter_plot_indiv.sh** (or **submit_filter_plot_indiv_regenie_UKBB.sh** for the UKBB GWAS, which includes UKBB-specific filtering steps) or in bulk across all studies with **batch_submit_filter_plot_indiv.R**. 

Scripts in *filter_public_sumstats/*:

The scripts in this folder are used to perform QC for the publicly downloaded summary statistics for some hormones.

- **add_allele2_from_ensembl.sh** - For summary statistics that only provide the effect allele, add in the other allele using VCFs from the appropriate Ensembl release by matching on chromosome, position, and effect allele.
- **add_MAF_from_1000G.sh** - For summary statistics that do not provide minor allele frequency (MAF), add in from the EUR-ancestry MAF in 1000 Genomes project VCFs by matching on chromosome, position, and both alleles.

The remaining scripts are similar to those in *filter_individual_results/* above, as they perform filtering and QC after the appropriate columns have been added and format summary statistics for meta-analysis in METAL.

Scripts in this folder to perform meta-analyses:

Example set of instructions for meta-analysis are in *meta_analysis_instructions/*, which can be submitted with the script **meta_all_shell.sh**. Variations on these (for example, meta-analysis excluding publicly available summary statistics, or meta-analysis excluding a single cohort) can be produced by changing the instructions in *meta_analysis_instructions/*.

After meta-analysis, **filter_plot_meta_results.R** (batch submitted across all strata with **batch_submit_filter_plot_meta_results.R**) is used to check the summary statistics for implausibly large standard errors and duplicate entries, as well as to plot Manhattan plots and QQ plots (split by MAF bin).

Scripts in *get_lead_snps/* to identify index variants after meta-analysis:

**get_lead_snps_distance_pruning.sh** - extracts genome-wide significant SNPs and windows of +/-500kb around these; also adds RSIDs to each genome-wide significant SNP. This script then calls **get_lead_snps_distance_pruning.R** to perform the distance-based pruning, i.e. to retain only the SNP with the lowest p-value within each genomic window. Pruning is done in two groups: (1) all SNPs, and (2) just SNPs in the 1000 Genomes EUR-ancestry subset, as this will be used to inform LD-based classification in subsequent steps. These jobs can be submitted in parallel across strata with **batch_submit_get_lead_snps.R**.

Scripts in *sensitivity_compare_public/* to compare lead variants from meta-analyses that include and don't include summary statistics from publicly available datasets.

