#!/bin/bash

cd /well/lindgren/samvida/hormones_infertility/exome_seq_results/conditioned_on_common

echo -e "phenotype\trare_variant\tcommon_snps" > snplist_to_condition.txt

declare -a HORMONES=("Oestradiol_F" "Testosterone_F" "Testosterone_M")

for hr in "${HORMONES[@]}"; do
	# Grab the lead variants on same chromosome and within +/-500 kb of each position
	
	head -n 1 /well/lindgren/samvida/hormones_infertility/meta_results_230613/lead_snps/all_lead_snp_sumstats_${hr}_EUR_with_rsids.txt \
	> nearby_common_variants/${hr}_nearby_common_variant_sumstats.txt

	tail -n +2 "sig_rare_variants/${hr}_signif_variants.txt" | while read -r variant_line; do
		chr_match=$(echo "$variant_line" | awk '{print $1}')
		# Change pattern_column1 from X to 23
		if [[ "$chr_match" == "X" ]]; then
			chr_match="23"
		fi
		pos_match=$(echo "$variant_line" | awk '{print $2}')
		
		awk -v pos_match="$pos_match" -v chr_match="$chr_match" '
			$3 == chr_match && $4 >= pos_match - 500000 && $4 <= pos_match + 500000
		' "/well/lindgren/samvida/hormones_infertility/meta_results_230613/lead_snps/all_lead_snp_sumstats_${hr}_EUR_with_rsids.txt" \
		> nearby_common_variants/tmp.txt
		
		cat nearby_common_variants/tmp.txt >> nearby_common_variants/${hr}_nearby_common_variant_sumstats.txt

		# Print SNP list per rare variant 
		rare_var=$(echo "$variant_line" | awk '{print $3}')
		common_snps=$(awk '{print $2}' nearby_common_variants/tmp.txt | paste -sd, -)
		echo -e "${hr}\t${rare_var}\t${common_snps}" >> snplist_to_condition.txt

		rm nearby_common_variants/tmp.txt
	done  
done

# Concatenate all the unique SNPs into a single file for plink
for fname in nearby_common_variants/*_nearby_common_variant_sumstats.txt; do
	# Extract the RSID, chromosome, and hg38 position
	tail -n +2 "$fname" | awk '{print $2 "\t" $3 "\t" $4}' \
	>> nearby_common_variants/all_rsids_to_get.txt
done
# Only retain the unique SNPs
sort nearby_common_variants/all_rsids_to_get.txt | uniq > nearby_common_variants/unique_rsids_to_get.txt
rm nearby_common_variants/all_rsids_to_get.txt
