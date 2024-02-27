Scripts to classify the identified hormone lead variants as reported or novel based on previously reported associations in the GWAS Catalog. Scripts must be run in the following order:

**get_gwascatalog_hormones_in_1kg_hg38.sh** is a script to subset the variants in GWAS Catalog associated with any of 28 reproductive hormones (list in **all_hormone_traits.txt**). It calls **make_gwascat_bed.R** to create a bed file for variant chromosome and position from the table of associations and then generates a list of all the 1000 Genomes variant ids so they can be compared to the SNPs from our study.

Scripts in *reported_novel_common_variants/*:

Steps for conditional analysis:

1. **make_sumstats_for_gcta.sh** - This needs to be run first to create files that are friendly for GCTA-COJO. Submit across all strata with **submit_make_sumstats_for_gcta.R**
2. **gcta_cojo_conditional_analysis.sh** - This is called by the script below that sorts variants into reported and unreported variants. For each unreported lead variant, GCTA-COJO first checks if all the SNPs in the window are collinear (which means the effect of the lead variant can be explained by GWAS Catalog variants), if not, it performs conditional analysis; if the lead variant is the only remaining variant, then it is novel, and if not we can check the conditional P-value to determine whether it is novel (P<0.05 even after conditioning on GWAS Catalog variants).

Scripts must be run in the following order on the output from lead SNP identification from GWAS meta-analyses.

1. **1a_sort_1kg_variants_for_conditional_analysis.R** - Decide which lead variants (from the 1000 Genomes project, as this is an LD-based analysis) need to be run through conditional analysis - any variant that has already been reported in the GWAS Catalog for this hormone does not need to be run - and create the 500 kb window of variants with previously reported hormone associations around each lead variant. Submit the script to perform conditional analysis, i.e. **gcta_cojo_conditional_analysis.sh**.
**1b_classify_1kg_variants.R** - Once GCTA-COJO has been run, check the "possibly novel" SNPs, i.e. check the conditional P-value to determine whether it is novel (P<0.05 even after conditioning on GWAS Catalog variants).
2. **2_classify_other_variants.R** - If a lead variant was not in the 1000 Genomes dataset or was on the X chromosome, we cannot use GCTA-COJO conditional analysis as we don't have LD information. Instead, we classify just based on distance: if a GWAS Catalog SNP is within 500kb of the lead variant, the lead variant is reported. 
3. **3_plots_and_tables.R** - Plot the number of reported/novel SNPs in each strata, and add functional annotations from the OpenTargets Catalog (such as nearest gene and v2g score) using **add_annotations_opentargets.sh**. 

Scripts in *reported_novel_rare_variants/*:

1. **1_sort_variants.R** - Across all rare variant and gene-based tests, check whether the associated gene has previously been reported for a hormone trait in GWAS Catalog. Plot the number of reported and novel genes in each strata.
