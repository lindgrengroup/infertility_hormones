#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 02/10/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J prep_betascan2
#SBATCH -o /well/lindgren/samvida/hormones_infertility/selection_analyses/logs/prep_betascan2-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/selection_analyses

mkdir data/tmp_betascan_sorted
# Subset BetaScan2 scores by chromosome and sort by position so it is ready for joining with sumstats
for CHR in {1..22}; do
	awk -v CHR="chr$CHR" -v OFS="\t" '{ if ($1 == CHR) { print $0 }}' \
	/well/lindgren/samvida/Resources/selection_scores_public/BetaScan2/hg38/all_B2std_hg38.txt \
	| sort -k 2 > data/tmp_betascan_sorted/chr${CHR}_sorted_beta_scores.txt
done

# Merge sumstats for infertility and betascan2 by chr, pos (hg38)
while IFS=$'\t' read -r -a STRATA_INFO
do
	for CHR in {1..22}; do
		# Infertility file columns: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, P.value, HetPVal, N_CASES, N_CONTROLS, maf, chr, pos
		# Hormone file columns: ID, CHROM, GENPOS, MAF, Allele1, Allele2, Freq1, FreqSE, BETA, SE, PVALUE, Direction, HetPVal
		awk -v CHR="$CHR" '{ if ($12 == CHR) { print $0 }}' \
		${STRATA_INFO[1]} | sort -k 13 > data/tmp_betascan_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt

		# Output file columns: CHR, POS, BetaScan2, BetaScan2_Std,
		# Tested_Allele_Freq, Tested_Allele, Other_Allele, BETA, PVALUE
		join -1 2 -2 13 -o 1.1,1.2,1.3,1.4,2.4,2.2,2.3,2.5,2.7 \
		<(sort -k 2 data/tmp_betascan_sorted/chr${CHR}_sorted_beta_scores.txt) \
		<(sort -k 13 data/tmp_betascan_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt) \
		> data/tmp_betascan_sorted/${STRATA_INFO[0]}_chr${CHR}_beta_scores.txt

		rm data/tmp_betascan_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt
	done
	cat data/tmp_betascan_sorted/${STRATA_INFO[0]}_chr*_beta_scores.txt >> data/${STRATA_INFO[0]}_betascan2_scores.txt
	# Add header column
	sed -i '1s/^/CHR POS BetaScan2 BetaScan2_Std Tested_Allele_Freq Tested_Allele Other_Allele BETA PVALUE\n/' \
	data/${STRATA_INFO[0]}_betascan2_scores.txt
done < infertility_files_list.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

