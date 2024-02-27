#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 2
#SBATCH -J add_allele2_from_ensembl
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/suhre_add_allele2_from_ensembl-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_ARRAY_TASK_ID" 
echo "Slurm Task ID: $SLURM_ARRAY_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/public_sumstats/tmp_Suhre_qc

for CHR in {1..22}; do
	# Create temporary vcf without the header and sort by position
	zgrep -v '^#' /well/lindgren/samvida/Resources/hg19/vcfs/homo_sapiens-chr${CHR}.vcf.gz | \
	sort -t $'\t' -k 2 -b > tmp_chr${CHR}_vcf_info.txt

	# Input files: 
	# (1) /well/lindgren/samvida/hormones_infertility/public_sumstats/FSH/filtered/Suhre_FSH_sex_comb_filtered.txt
	# (2) /well/lindgren/samvida/hormones_infertility/public_sumstats/LH/filtered/Suhre_LH_sex_comb_filtered.txt

	awk -v CHR="$CHR" '{ if ($2==CHR) { print } }' \
	/well/lindgren/samvida/hormones_infertility/public_sumstats/FSH/filtered/Suhre_FSH_sex_comb_filtered.txt | \
	sort -t $'\t' -k 3 -b > tmp_FSH_chr${CHR}_sumstats.txt
	awk -v CHR="$CHR" '{ if ($2==CHR) { print } }' \
	/well/lindgren/samvida/hormones_infertility/public_sumstats/LH/filtered/Suhre_LH_sex_comb_filtered.txt | \
	sort -t $'\t' -k 3 -b > tmp_LH_chr${CHR}_sumstats.txt

	# Merge in the ref/alt alleles from vcf by position
	join -t $'\t' -1 3 -2 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.4,2.5 \
	tmp_FSH_chr${CHR}_sumstats.txt tmp_chr${CHR}_vcf_info.txt > tmp_FSH_chr${CHR}_check_alleles.txt

	join -t $'\t' -1 3 -2 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.4,2.5 \
	tmp_LH_chr${CHR}_sumstats.txt tmp_chr${CHR}_vcf_info.txt > tmp_LH_chr${CHR}_check_alleles.txt

	# Replace ALLELE0 (NA) with the non-matching allele
	awk '{ if ($4==$13) { $5=$14 } else if ($4==$14) { $5=$13 } else { $4=$5="NA" }} { print }' \
	tmp_FSH_chr${CHR}_check_alleles.txt > tmp_FSH_chr${CHR}_sumstats_merged.txt
	awk '{ if ($4==$13) { $5=$14 } else if ($4==$14) { $5=$13 } else { $4=$5="NA" }} { print }' \
	tmp_LH_chr${CHR}_check_alleles.txt > tmp_LH_chr${CHR}_sumstats_merged.txt

	# Remove any rows where both alleles are now NA as that doesn't match the vcf
	awk '!($4=="NA" && $5=="NA")' tmp_FSH_chr${CHR}_sumstats_merged.txt > FSH_chr${CHR}_sumstats.txt
	awk '!($4=="NA" && $5=="NA")' tmp_LH_chr${CHR}_sumstats_merged.txt > LH_chr${CHR}_sumstats.txt

	rm tmp_*
done

# After all chromosomes have run
# Concatenate results across chromosomes
cd /well/lindgren/samvida/hormones_infertility/public_sumstats
head -n 1 FSH/filtered/og_Suhre_FSH_sex_comb_filtered.txt > FSH/filtered/Suhre_FSH_sex_comb_filtered.txt
head -n 1 LH/filtered/og_Suhre_LH_sex_comb_filtered.txt > LH/filtered/Suhre_LH_sex_comb_filtered.txt
for CHR in {1..22}; do
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' tmp_Suhre_qc/FSH_chr${CHR}_sumstats.txt >> FSH/filtered/Suhre_FSH_sex_comb_filtered.txt
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' tmp_Suhre_qc/LH_chr${CHR}_sumstats.txt >> LH/filtered/Suhre_LH_sex_comb_filtered.txt
done

rm -r tmp_Suhre_qc

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
