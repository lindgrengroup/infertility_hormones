#!/bin/bash

# Variant-level results
cd /well/lindgren/samvida/hormones_infertility/exome_seq_results/results_2310
gunzip *.gz

declare -a STRATA=("FSH_F" "FSH_M" "LH_F" "LH_M" "Oestradiol_F" "Oestradiol_M" "Testosterone_F" "Testosterone_M" "female_infertility_binary" "idiop_infertility_exclusion_binary" "male_infertility_binary")

for strata_name in "${STRATA[@]}"; do
	# Get variant-level results
	header=$(head -n 1 chr1_${strata_name}_EUR.txt.singleAssoc.txt)
	echo "${header}" > cat_results/${strata_name}_variant.txt
	tail -q -n +2 chr*_${strata_name}_EUR.txt.singleAssoc.txt >> cat_results/${strata_name}_variant.txt
	# Get gene-level results
	header=$(head -n 1 chr1_${strata_name}_EUR.txt)
	echo "${header}" > cat_results/${strata_name}_gene.txt
	tail -q -n +2 chr*_${strata_name}_EUR.txt >> cat_results/${strata_name}_gene.txt	
done
