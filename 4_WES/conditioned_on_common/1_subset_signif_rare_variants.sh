#!/bin/bash

cd /well/lindgren/samvida/hormones_infertility/exome_seq_results

declare -a STRATA=("FSH_F" "FSH_M" "LH_F" "LH_M" "Oestradiol_F" "Oestradiol_M" "Testosterone_F" "Testosterone_M" "female_infertility_binary" "idiop_infertility_exclusion_binary" "male_infertility_binary")

for strata_name in "${STRATA[@]}"; do
	header=$(head -n 1 results_2310/cat_results/${strata_name}_variant.txt)
	echo "${header}" > conditioned_on_common/sig_rare_variants/${strata_name}_signif_variants.txt
	awk '$13 < 1E-07' results_2310/cat_results/${strata_name}_variant.txt >> \
	conditioned_on_common/sig_rare_variants/${strata_name}_signif_variants.txt

	# Delete the file if there were no significant rare variants
	[ $(wc -l < conditioned_on_common/sig_rare_variants/${strata_name}_signif_variants.txt) -eq 1 ] \
	&& rm conditioned_on_common/sig_rare_variants/${strata_name}_signif_variants.txt
done

