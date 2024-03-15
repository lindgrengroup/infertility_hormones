Scripts to test for directional selection as measured by CMS, aDNA, and SDS at infertility or hormone loci and vice-versa; also to test for balancing selection using BetaScan2 scores. Protocol first described by Mathieson et al. 2023 - https://pubmed.ncbi.nlm.nih.gov/36864135/.

For aDNA (Mathieson et al. 2015) and CMS (Grossman et al. 2013), we don't have access to variant level scores so we just check what the minimum P-values are for infertility/hormone traits in loci identified to be under selection. Do this using **0_get_sumstats_in_selection_loci.sh** to extract summary statistics at the loci of interest, and plot and table the results with **0_table_plot_infertility_sumstats_in_selection_loci.R**.

For BetaScan2 (Siewert et al. 2020) and SDS (Field et al. 2016), we have variant level scores, so we can compare directly to the infertility/hormone summary statistics.

1. **1_prepare_selection_scores.R** - SDS has to be converted to a P-value; BetaScan2 scores have to be normalised and create a P-value.
2. **2_liftover_selection_scores_hg38.sh** - Selection scores are on the hg19 genome build, lift over to hg38.
3. **3_table_plot_selection_scores_in_infertility_loci.R** - Assess minimum/maximum selection scores around windows defined by infertility lead variants.
4. **4_prep_trait_sds.sh** - For analyses where SDS is aligned to the trait-increasing allele rather than the derived allele, so we can test whether the direction of selection is for or against the trait. Plot these using **4_trait_sds_scores_plots.R**.
	**4_prep_betascan2_scores.sh** - For analyses where the BetaScan2 score is added to the GWAS summary statistics so we can see whether the same variants have high scores for both GWAS and balancing selection. Plot these using **4_betascan2_plots.R**.
