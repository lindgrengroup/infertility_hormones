#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 18/08/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J sensitivity_age_menopause
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/sensitivity_age_menopause-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS
mkdir sensitivity_age_menopause

declare -a STRATA=("FSH" "LH")

for hormone in "${STRATA[@]}"; do
	for CHR in {1..22}; do
		# Sort the files by rsid to join
		sort -k 3 REGENIE_results/${hormone}_non_finnish_EUR_F/step2/step2_out_${hormone}_non_finnish_EUR_F_chr${CHR}_${hormone}.regenie \
		-o REGENIE_results/${hormone}_non_finnish_EUR_F/step2/step2_out_${hormone}_non_finnish_EUR_F_chr${CHR}_${hormone}.regenie
		sort -k 3 REGENIE_results/${hormone}_normal_menopause/step2/step2_out_${hormone}_normal_menopause_chr${CHR}_${hormone}.regenie \
		-o REGENIE_results/${hormone}_normal_menopause/step2/step2_out_${hormone}_normal_menopause_chr${CHR}_${hormone}.regenie
		# Join both sets of results
		join -1 3 -2 3 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.8,1.10,1.11,1.13,2.8,2.10,2.11,2.13 \
		REGENIE_results/${hormone}_non_finnish_EUR_F/step2/step2_out_${hormone}_non_finnish_EUR_F_chr${CHR}_${hormone}.regenie \
		REGENIE_results/${hormone}_normal_menopause/step2/step2_out_${hormone}_normal_menopause_chr${CHR}_${hormone}.regenie \
		> sensitivity_age_menopause/${hormone}_chr${CHR}_joined.txt

		# Subset to only results that are sub-significant P<1E-05 (-log10P>5) in either GWAS
		awk '($10 ~ /^[0-9]+(\.[0-9]+)?$/ && $14 ~ /^[0-9]+(\.[0-9]+)?$/) && ($10 > 5 || $14 > 5)' sensitivity_age_menopause/${hormone}_chr${CHR}_joined.txt \
		> sensitivity_age_menopause/signif_${hormone}_chr${CHR}_joined.txt

        # Subset to only lead SNPs from the full meta-analyses
        awk 'NR > 1 {print $2}' /well/lindgren/samvida/hormones_infertility/meta_results_230613/lead_snps/all_lead_snp_sumstats_${hormone}_F_all_anc_with_rsids.txt \
        > sensitivity_age_menopause/grep_${hormone}_rsids.txt
        grep -Fwf sensitivity_age_menopause/grep_${hormone}_rsids.txt sensitivity_age_menopause/signif_${hormone}_chr${CHR}_joined.txt \
        > sensitivity_age_menopause/MA_lead_snps_${hormone}_chr${CHR}.txt
	done

    # ChrX
    CHR="X"
    # Sort the files by rsid to join
	sort -k 3 REGENIE_results/${hormone}_non_finnish_EUR_F/step2/step2_out_${hormone}_non_finnish_EUR_F_chr${CHR}_${hormone}.regenie \
	-o REGENIE_results/${hormone}_non_finnish_EUR_F/step2/step2_out_${hormone}_non_finnish_EUR_F_chr${CHR}_${hormone}.regenie
	sort -k 3 REGENIE_results/${hormone}_normal_menopause/step2/step2_out_${hormone}_normal_menopause_chr${CHR}_${hormone}.regenie \
	-o REGENIE_results/${hormone}_normal_menopause/step2/step2_out_${hormone}_normal_menopause_chr${CHR}_${hormone}.regenie
	# Join both sets of results
	join -1 3 -2 3 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.8,1.10,1.11,1.13,2.8,2.10,2.11,2.13 \
	REGENIE_results/${hormone}_non_finnish_EUR_F/step2/step2_out_${hormone}_non_finnish_EUR_F_chr${CHR}_${hormone}.regenie \
	REGENIE_results/${hormone}_normal_menopause/step2/step2_out_${hormone}_normal_menopause_chr${CHR}_${hormone}.regenie \
	> sensitivity_age_menopause/${hormone}_chr${CHR}_joined.txt

	# Subset to only results that are sub-significant P<1E-05 (-log10P>5) in either GWAS
	awk '($10 ~ /^[0-9]+(\.[0-9]+)?$/ && $14 ~ /^[0-9]+(\.[0-9]+)?$/) && ($10 > 5 || $14 > 5)' sensitivity_age_menopause/${hormone}_chr${CHR}_joined.txt \
	> sensitivity_age_menopause/signif_${hormone}_chr${CHR}_joined.txt

    # Subset to only lead SNPs from the full meta-analyses
    awk 'NR > 1 {print $2}' /well/lindgren/samvida/hormones_infertility/meta_results_230613/lead_snps/all_lead_snp_sumstats_${hormone}_F_all_anc_with_rsids.txt \
    > sensitivity_age_menopause/grep_${hormone}_rsids.txt
    grep -Fwf sensitivity_age_menopause/grep_${hormone}_rsids.txt sensitivity_age_menopause/signif_${hormone}_chr${CHR}_joined.txt \
    > sensitivity_age_menopause/MA_lead_snps_${hormone}_chr${CHR}.txt

	# Concatenate across chromosomes with header
	echo -e "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N_all BETA_all SE_all LOG10P_all N_meno BETA_meno SE_meno LOG10P_meno" \
    > sensitivity_age_menopause/signif_${hormone}.txt
	cat sensitivity_age_menopause/signif_${hormone}_chr*_joined.txt >> sensitivity_age_menopause/signif_${hormone}.txt
    
    # Concatenate across chromosomes with header
	echo -e "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N_all BETA_all SE_all LOG10P_all N_meno BETA_meno SE_meno LOG10P_meno" \
    > sensitivity_age_menopause/MA_lead_snps_${hormone}.txt
	cat sensitivity_age_menopause/MA_lead_snps_${hormone}_chr*.txt >> sensitivity_age_menopause/MA_lead_snps_${hormone}.txt
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
