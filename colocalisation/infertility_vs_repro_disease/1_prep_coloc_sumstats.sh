#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J prep_coloc_sumstats
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/prep_coloc_sumstats-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/colocalisation

# Lift reproductive trait sumstats up to hg38 and add sample sizes
module load liftOver/20210519

while IFS=$'\t' read -r -a STRATA_INFO
do
	# Replace 23 with X for chrX
	awk -v FS="\t" -v OFS="\t" '{ sub("23", "X", $2) } 1' \
	/well/lindgren-ukbb/projects/ukbb-11867/samvida/obesity_wrh_thesis/two_sample_mr/repro_outcomes/${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	> sentinel_windows/liftover/tmp_${STRATA_INFO[0]}.txt
	
	# Create bed file with chrX, pos0, pos1, RSID
	awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$2, $3-1, $3, $1}' \
	sentinel_windows/liftover/tmp_${STRATA_INFO[0]}.txt \
	> sentinel_windows/liftover/${STRATA_INFO[0]}_hg19.bed
			
	liftOver sentinel_windows/liftover/${STRATA_INFO[0]}_hg19.bed \
	/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
	sentinel_windows/liftover/${STRATA_INFO[0]}_hg38.bed \
	sentinel_windows/liftover/${STRATA_INFO[0]}_unlifted.bed

	# Create a new sumstats file by joining via the RSID and add sample sizes for Ncas and Nctrl
	# unique_id: 4th column in the lifted over bed file and 1st column in the original sumstats
	join -1 4 -2 1 -o 2.1,1.1,1.3,2.4,2.5,2.6,2.7,2.8,2.9 \
	<(sort -k 4 sentinel_windows/liftover/${STRATA_INFO[0]}_hg38.bed) \
	<(sort -k 1 /well/lindgren-ukbb/projects/ukbb-11867/samvida/obesity_wrh_thesis/two_sample_mr/repro_outcomes/${STRATA_INFO[0]}_sumstats_with_rsids.txt) \
	> sentinel_windows/liftover/${STRATA_INFO[0]}_hg38.txt
	awk -v NCAS="${STRATA_INFO[1]}" -v NCON="${STRATA_INFO[2]}" 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, NCAS, NCON }' \
	sentinel_windows/liftover/${STRATA_INFO[0]}_hg38.txt \
	> sentinel_windows/liftover/${STRATA_INFO[0]}_hg38_with_samplesizes.txt

	sed -i '1i\RSID\tCHR\tPOS\tALLELE1\tALLELE2\tFREQ1\tBETA\tSE\tPVALUE\tNCASE\tNCONTROL' \
	sentinel_windows/liftover/${STRATA_INFO[0]}_hg38_with_samplesizes.txt

done < repro_files_list.txt

rm sentinel_windows/liftover/tmp_*
rm sentinel_windows/liftover/*.bed

# Get windows of sumstats (+/- 50kb) around each infertility lead SNP
for finf in female_infertility_analysis1 female_infertility_analysis3 female_infertility_analysis4 female_infertility_analysis5; do
	# Read lead SNPs
	while IFS=$'\t' read -r -a SNP_INFO; do
		awk -v CHR="${SNP_INFO[1]}" -v POS="${SNP_INFO[2]}" '$12==CHR && $13>=(POS-50000) && $13<=(POS+50000)' \
		/well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/${finf}_eur_MA_results_chr_pos.txt \
		> sentinel_windows/${SNP_INFO[0]}_${finf}.txt
		sed -i '1i\MarkerName\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tP.value\tHetPVal\tN_CASES\tN_CONTROLS\tmaf\tchr\tpos' \
		sentinel_windows/${SNP_INFO[0]}_${finf}.txt

		for repro_trait in endometriosis hmb PCOS uterine_fibroids; do
			awk -v CHR="${SNP_INFO[1]}" -v POS="${SNP_INFO[2]}" '$2=="chr"CHR && $3>=(POS-50000) && $3<=(POS+50000)' \
			sentinel_windows/liftover/${repro_trait}_hg38_with_samplesizes.txt \
			> sentinel_windows/${SNP_INFO[0]}_${repro_trait}.txt
			sed -i '1i\RSID\tCHR\tPOS\tALLELE1\tALLELE2\tFREQ1\tBETA\tSE\tPVALUE\tNCASE\tNCONTROL' \
			sentinel_windows/${SNP_INFO[0]}_${repro_trait}.txt
		done
	done < infertility_sentinels/${finf}_eur_sentinels.txt
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
