Scripts to assess local genetic correlations between: infertility, reproductive hormones, reproductive diseases, and obesity. Only run traits with significant heritability Z>4 as we need accurate gcov_int estimates. Scripts adapted from https://github.com/josefin-werme/LAVA. 

Pre-processing:

**0_extract_cross_trait_rg_ldsc.sh** - Uses Duncan's modified LDSC https://github.com/astheeggeggs/UKBB_ldsc_r2 to get pairwise genetic correlations among all traits of interest, because we want to feed the sample overlap (gcov_int) into LAVA. Also prepare the locus file.
**0_summarise_cross_trait_ldsc.R** - Create matrix of correlations from the cross-trait gcov_int column.
**1_prepare_blocks_loci.sh** - Split the blocks file from LAVA developers into a separate list of loci per chromosome.

Analysis: 
**2_run_lava_bivariate.R** - For a given pair of phenotypes, check univariate heritability at each locus and run bivariate rG if significant. Parallelise this across chromosomes and run across every pairwise combination of traits using **2_submit_lava_bivariate.sh** and **2_batch_submit_lava_bivariate.R** 
**3_plot_bivariate_results.R** - Manhattan-like plot across all phenotype-pairs tested; counts plot to summarise the number of locally significant rG per pair.
**4_run_lava_multivariate.R** - At loci with multiple associated phenotypes, run multiple regression.