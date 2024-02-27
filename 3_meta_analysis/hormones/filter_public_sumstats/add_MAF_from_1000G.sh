#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J add_MAF_from_1kg
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/suhre_add_AF_from_1kg-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_ARRAY_TASK_ID" 
echo "Slurm Task ID: $SLURM_ARRAY_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility

mkdir tmp_Suhre_Prins_qc

for CHR in {1..23}; do
	if [[ "$CHR" -eq 23 ]]; then
		chrvcf="X"
	else
		chrvcf=$CHR
	fi

	cd /well/lindgren/samvida/hormones_infertility

	# Create temporary vcf without the header and sort by position
	zgrep -v '^#' /well/lindgren/samvida/Resources/hg19/vcfs/1000GENOMES-phase_3.vcf.gz | \
	awk -v chrvcf="$chrvcf" '{ if ($1==chrvcf) { print } }' - \
	> sorted_chr_vcfs_for_munging/tmp_chr${chrvcf}_1kg_vcf_info.txt

	sort -t $'\t' -k 2 -b sorted_chr_vcfs_for_munging/tmp_chr${chrvcf}_1kg_vcf_info.txt \
	> sorted_chr_vcfs_for_munging/chr${chrvcf}_1kg_vcf_info.txt
	
	# Retain MAF for EUR, which is specified as ";EUR=...;" 
	awk -v "OFS=\t" '{$8=$8;sub(/.*EUR=/, "", $8); sub(/;.*/, "", $8); print}' \
	sorted_chr_vcfs_for_munging/chr${chrvcf}_1kg_vcf_info.txt \
	> sorted_chr_vcfs_for_munging/chr${chrvcf}_1kg_vcf_EUR_MAF.txt

	# SUHRE ET AL.
	# Input files: 
	# (1) /well/lindgren/samvida/hormones_infertility/public_sumstats/FSH/filtered/Suhre_FSH_sex_comb_filtered.txt
	# (2) /well/lindgren/samvida/hormones_infertility/public_sumstats/LH/filtered/Suhre_LH_sex_comb_filtered.txt
	# PRINS ET AL.
	# (3) /well/lindgren/samvida/hormones_infertility/public_sumstats/Testosterone/filtered/og_Prins_Testosterone_M_filtered.txt

	awk -v CHR="$CHR" '{ if ($2==CHR) { print } }' \
	./public_sumstats/FSH/filtered/og_Suhre_FSH_sex_comb_filtered.txt | \
	sort -t $'\t' -k 3 -b > tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_sumstats.txt
	
	awk -v CHR="$CHR" '{ if ($2==CHR) { print } }' \
	./public_sumstats/LH/filtered/og_Suhre_LH_sex_comb_filtered.txt | \
	sort -t $'\t' -k 3 -b > tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_sumstats.txt
	
	awk -v CHR="$CHR" '{ if ($2==CHR) { print } }' \
	./public_sumstats/Testosterone/filtered/og_Prins_Testosterone_M_filtered.txt | \
	sort -t $'\t' -k 3 -b > tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_sumstats.txt

	# Merge in the ref/alt alleles and MAF from vcf by position
	join -t $'\t' -1 3 -2 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.8,1.8,1.9,1.10,1.11,1.12,2.4,2.5 \
	tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_sumstats.txt sorted_chr_vcfs_for_munging/chr${chrvcf}_1kg_vcf_EUR_MAF.txt \
	> tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_check_alleles.txt

	join -t $'\t' -1 3 -2 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.8,1.8,1.9,1.10,1.11,1.12,2.4,2.5 \
	tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_sumstats.txt sorted_chr_vcfs_for_munging/chr${chrvcf}_1kg_vcf_EUR_MAF.txt \
	> tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_check_alleles.txt
	
	join -t $'\t' -1 3 -2 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.8,1.8,1.9,1.10,1.11,1.12,2.4,2.5 \
	tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_sumstats.txt sorted_chr_vcfs_for_munging/chr${chrvcf}_1kg_vcf_EUR_MAF.txt \
	> tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_check_alleles.txt

	# Only retain rows with A/G/C/T ref and alt alleles
	awk -v "OFS=\t" '{ if ($13 ~ "^A$|^C$|^G$|^T$" && $14 ~ "^A$|^C$|^G$|^T$") { print }}' \
	tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_check_alleles.txt > tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_acgt.txt
	awk -v "OFS=\t" '{ if ($13 ~ "^A$|^C$|^G$|^T$" && $14 ~ "^A$|^C$|^G$|^T$") { print }}' \
	tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_check_alleles.txt > tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_acgt.txt
	awk -v "OFS=\t" '{ if ($13 ~ "^A$|^C$|^G$|^T$" && $14 ~ "^A$|^C$|^G$|^T$") { print }}' \
	tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_check_alleles.txt > tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_acgt.txt

	# Column order: $4 is ALLELE1, $5 is ALLELE0, $13 is REF, and $14 is ALT
	# $6 is A1FREQ, $7 is MAF

	# For Suhre et al., Replace ALLELE0 (NA) with the non-matching allele
	# And A1FREQ with MAF or 1-MAF depending on whether ALLELE1 matches ALT or not
	awk -v "OFS=\t" '{ if ($4==$13) { $5=$14; $6=1-$7 } else if ($4==$14) { $5=$13; $6=$7 } else { $4=$5=$6="NA" }} { print }' \
	tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_acgt.txt > tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_sumstats_merged.txt
	
	awk -v "OFS=\t" '{ if ($4==$13) { $5=$14; $6=1-$7 } else if ($4==$14) { $5=$13; $6=$7 } else { $4=$5=$6="NA" }} { print }' \
	tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_acgt.txt > tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_sumstats_merged.txt

	# For Prins et al., check that both ALLELE0 and ALLELE1 match the ref/alt 
	# And replace A1FREQ with MAF or 1-MAF depending on whether ALLELE1 matches ALT or not
	awk -v "OFS=\t" '{ if ($4==$13 && $5==$14) { $6=1-$7 } else if ($4==$14 && $5==$13) { $6=$7 } else { $4=$5=$6="NA" }} { print }' \
	tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_acgt.txt > tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_sumstats_merged.txt
	
	# Remove any rows where both alleles are now NA as that doesn't match the vcf
	awk '!($4=="NA" && $5=="NA")' tmp_Suhre_Prins_qc/tmp_FSH_chr${CHR}_sumstats_merged.txt > tmp_Suhre_Prins_qc/FSH_chr${CHR}_sumstats.txt
	awk '!($4=="NA" && $5=="NA")' tmp_Suhre_Prins_qc/tmp_LH_chr${CHR}_sumstats_merged.txt > tmp_Suhre_Prins_qc/LH_chr${CHR}_sumstats.txt
	awk '!($4=="NA" && $5=="NA")' tmp_Suhre_Prins_qc/tmp_Testosterone_chr${CHR}_sumstats_merged.txt > tmp_Suhre_Prins_qc/Testosterone_chr${CHR}_sumstats.txt

	rm tmp_Suhre_Prins_qc/tmp_*
done

# After all chromosomes have run
# Concatenate results across chromosomes
head -n 1 public_sumstats/FSH/filtered/og_Suhre_FSH_sex_comb_filtered.txt > public_sumstats/FSH/filtered/Suhre_FSH_sex_comb_filtered.txt
head -n 1 public_sumstats/LH/filtered/og_Suhre_LH_sex_comb_filtered.txt > public_sumstats/LH/filtered/Suhre_LH_sex_comb_filtered.txt
head -n 1 public_sumstats/Testosterone/filtered/og_Prins_Testosterone_M_filtered.txt > public_sumstats/Testosterone/filtered/Prins_Testosterone_M_filtered.txt

for CHR in {1..23}; do
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' tmp_Suhre_Prins_qc/FSH_chr${CHR}_sumstats.txt >> public_sumstats/FSH/filtered/Suhre_FSH_sex_comb_filtered.txt
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' tmp_Suhre_Prins_qc/LH_chr${CHR}_sumstats.txt >> public_sumstats/LH/filtered/Suhre_LH_sex_comb_filtered.txt
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' tmp_Suhre_Prins_qc/Testosterone_chr${CHR}_sumstats.txt >> public_sumstats/Testosterone/filtered/Prins_Testosterone_M_filtered.txt
done

rm -r tmp_Suhre_Prins_qc

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
