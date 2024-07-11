Scripts to assess local genetic correlations between: infertility, reproductive hormones, reproductive diseases, and obesity. 

Pre-processing:

**0_extract_cross_trait_rg_ldsc.sh** - Uses Duncan's modified LDSC https://github.com/astheeggeggs/UKBB_ldsc_r2 to get pairwise genetic correlations among all traits of interest, because we want to feed the sample overlap (gcov_int) into LAVA.