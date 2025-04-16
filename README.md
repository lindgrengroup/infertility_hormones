# Genome-wide analyses identify 25 infertility loci and relationships with reproductive traits across the allele frequency spectrum

Code for analyses in published article here - https://www.nature.com/articles/s41588-025-02156-8

Cite as: Venkatesh, S.S., Wittemans, L.B.L., Palmer, D.S. et al. Genome-wide analyses identify 25 infertility loci and relationships with reproductive traits across the allele frequency spectrum. Nat Genet (2025). https://doi.org/10.1038/s41588-025-02156-8

## Structure

- *1_phenotype_identification/* - Identify and perform quality control of binary traits (infertility cases and controls) and quantitative traits (hormone values) in the UK Biobank. 
- *2_GWAS/* - Perform genome-wide association studies across autosomes and the X chromosome, for binary traits and quantitative traits in the UK Biobank, using REGENIE. 
- *3_meta_analysis/* - Filter and QC individual cohort summary statistics, perform meta-analyses, filter, QC and plot results, and get lead SNPs. 
- *4_WES/* - Plot rare-variant and gene-based test results from whole exome sequencing (WES) analyses in the UK Biobank. Scripts to perform WES testing is all documented in Duncan Palmer's repository here - https://github.com/astheeggeggs/BRaVa_curation.
- *characterise_loci/* - Classify the identified common and rare variants (for hormones) as reported or novel based on previously published associations in the GWAS Catalog. 
- *colocalisation* - Check for statistical colocalisation between pairs of traits (such as infertility and hormones, or infertility and reproductive disease) as well as between our traits and GTeX eQTLs to identify putative causal genes.
- *combine_wes_gwas* - Compare effect sizes from GWAS meta-analyses with UK Biobank WES analyses for genes reported in both sets of analyses. 
- *heritability_rg* - LDSC-based analyses to calculate trait heritability, genetic correlations between pairs of traits (including reproductive diseases and the Neale lab UKB-wide phenotypes), and partitioned heritability analyses for tissue-specific and cell-type-specific enrichment. 
- *lava_local_rg* - LAVA analyses to calculate local genetic correlation between pairs of traits; added in response to reviewer request.
- *mixer_overlap* - MiXeR analyses to calculate polygenic overlap, regardless of genetic correlation, between pairs of traits; added in response to reviewer request.
- *selection_analyses* - Test for evolutionary selection as measured by CMS, aDNA, and SDS at infertility or hormone loci and vice-versa.
- *two_sample_mr* - Mendelian randomisation to test for genetically predicted causal effects of obesity and hormone levels on infertility.
