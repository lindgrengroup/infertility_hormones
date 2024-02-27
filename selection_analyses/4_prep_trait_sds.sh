#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 02/10/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J prep_trait_sds
#SBATCH -o /well/lindgren/samvida/hormones_infertility/selection_analyses/logs/prep_trait_sds-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/selection_analyses

mkdir data/tmp_sds_sorted
# Subset SDS scores by chromosome and sort by position so it is ready for joining with sumstats
for CHR in {1..22}; do
	awk -v CHR="chr$CHR" -v OFS="\t" '{ if ($1 == CHR) { print $0 }}' \
	/well/lindgren/samvida/Resources/selection_scores_public/SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals_hg38.txt \
	| sort -k 2 > data/tmp_sds_sorted/chr${CHR}_sorted_sds_scores.txt
done

# Merge sumstats for hormones and SDS by chr, pos (hg38)
while IFS=$'\t' read -r -a STRATA_INFO
do
	for CHR in {1..22}; do
		# Hormone file columns: ID, CHROM, GENPOS, MAF, Allele1, Allele2, Freq1, FreqSE, BETA, SE, PVALUE, Direction, HetPVal
		awk -v CHR="$CHR" '{ if ($2 == CHR) { print $0 }}' \
		${STRATA_INFO[1]} | sort -k 3 > data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt

		# Output file columns: CHR, POS, RSID, Ancestral_Allele, Derived_Allele, Derived_Allele_Freq, SDS,
		# Tested_Allele_Freq, Tested_Allele, Other_Allele, BETA, PVALUE
		join -1 2 -2 3 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.7,2.5,2.6,2.9,2.11 \
		<(sort -k 2 data/tmp_sds_sorted/chr${CHR}_sorted_sds_scores.txt) \
		<(sort -k 3 data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt) \
		> data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt

		# Remove rows where the AA/DA don't match Tested/Other
		awk -v FS=" " -v OFS="\t" '{ if ($4 == $9 && $5 == $10 || $4 == $10 && $5 == $9) { print $0 } }' \
		data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt \
		> data/tmp_sds_sorted/${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt
		
		rm data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt
		rm data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt
	done
	cat data/tmp_sds_sorted/${STRATA_INFO[0]}_chr*_sds_scores.txt >> data/${STRATA_INFO[0]}_sds_scores.txt
	# Add header column
	sed -i '1s/^/CHR POS RSID Ancestral_Allele Derived_Allele Derived_Allele_Freq SDS Tested_Allele_Freq Tested_Allele Other_Allele BETA PVALUE\n/' \
	data/${STRATA_INFO[0]}_sds_scores.txt
done < hormone_files_list.txt

# Merge sumstats for hormones and SDS by chr, pos (hg38)
while IFS=$'\t' read -r -a STRATA_INFO
do
	for CHR in {1..22}; do
		# Infertility file columns: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, P.value, HetPVal, N_CASES, N_CONTROLS, maf, chr, pos
		# Hormone file columns: ID, CHROM, GENPOS, MAF, Allele1, Allele2, Freq1, FreqSE, BETA, SE, PVALUE, Direction, HetPVal
		awk -v CHR="$CHR" '{ if ($12 == CHR) { print $0 }}' \
		${STRATA_INFO[1]} | sort -k 13 > data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt

		# Output file columns: CHR, POS, RSID, Ancestral_Allele, Derived_Allele, Derived_Allele_Freq, SDS,
		# Tested_Allele_Freq, Tested_Allele, Other_Allele, BETA, PVALUE
		join -1 2 -2 13 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.4,2.2,2.3,2.5,2.7 \
		<(sort -k 2 data/tmp_sds_sorted/chr${CHR}_sorted_sds_scores.txt) \
		<(sort -k 13 data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt) \
		> data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt

		# Remove rows where the AA/DA don't match Tested/Other
		awk -v FS=" " -v OFS=" " '{ if ($4 == toupper($9) && $5 == toupper($10) || $4 == toupper($10) && $5 == toupper($9)) { print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9),toupper($10),$11,$12 } }' \
		data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt \
		> data/tmp_sds_sorted/${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt
		
		rm data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sorted.txt
		rm data/tmp_sds_sorted/tmp_${STRATA_INFO[0]}_chr${CHR}_sds_scores.txt
	done
	cat data/tmp_sds_sorted/${STRATA_INFO[0]}_chr*_sds_scores.txt >> data/${STRATA_INFO[0]}_sds_scores.txt
	# Add header column
	sed -i '1s/^/CHR POS RSID Ancestral_Allele Derived_Allele Derived_Allele_Freq SDS Tested_Allele_Freq Tested_Allele Other_Allele BETA PVALUE\n/' \
	data/${STRATA_INFO[0]}_sds_scores.txt
done < infertility_files_list.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

