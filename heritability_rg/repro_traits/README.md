Scripts to assess rG of traits in our study against TSH, AMH, and dizygotic twinning and female reproductive diseases - endometriosis, PCOS, heavy menstrual bleeding, and uterine fibroids. 

- **external_sumstats_munge_calc_rg.sh** - Prepare summary statistics for rG analyses using LDSC. External sumstats are for TSH (Williams et al. 2023), AMH (Verdiesen et al. 2022), and dizygotic twinning (Mbarek et al. 2024), which were downloaded from GWAS Catalog. Use LDSC to calculate rG between all pairs of traits. Only report rG for traits with significant heritability (Z>4). 
- **munge_sumstats_calc_rg.sh** - Prepare summary statistics for rG analyses using LDSC. Reproductive disease sumstats are from in-house meta-analyses, which will be released on GWAS Catalog. Use LDSC to calculate rG between all pairs of traits. Only report rG for traits with significant heritability (Z>4). 
- **plot_rg.R** - Create the plots in the manuscript.
