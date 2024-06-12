#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J prep_coloc_sumstats_age_at_menopause
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/prep_coloc_sumstats_age_at_menopause-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/colocalisation

# Lift reproductive trait sumstats up to hg38 and add sample sizes
module load liftOver/20210519

# Replace 23 with X for chrX
awk -v FS="\t" -v OFS="\t" '{ sub("23", "X", $1) } 1' \
/well/lindgren/samvida/hormones_infertility/age_at_menopause/GCST90292553.tsv \
> sentinel_windows/liftover/tmp_age_at_menopause.txt
	
# Create bed file with chrX, pos0, pos1, RSID
awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$1, $2-1, $2, $9}' \
sentinel_windows/liftover/tmp_age_at_menopause.txt \
> sentinel_windows/liftover/age_at_menopause_hg19.bed
			
liftOver sentinel_windows/liftover/age_at_menopause_hg19.bed \
/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
sentinel_windows/liftover/age_at_menopause_hg38.bed \
sentinel_windows/liftover/age_at_menopause_unlifted.bed

# Create a new sumstats file by joining via the RSID and add sample sizes for N
# unique_id: 4th column in the lifted over bed file and 9th column in the original sumstats
join -1 4 -2 9 -o 2.9,1.1,1.3,2.3,2.4,2.5,2.6,2.7,2.8 \
<(sort -k 4 sentinel_windows/liftover/age_at_menopause_hg38.bed) \
<(sort -k 9 /well/lindgren/samvida/hormones_infertility/age_at_menopause/GCST90292553.tsv) \
> sentinel_windows/liftover/age_at_menopause_hg38.txt
awk -v N="173424" 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, N }' \
sentinel_windows/liftover/age_at_menopause_hg38.txt \
> sentinel_windows/liftover/age_at_menopause_hg38_with_samplesizes.txt

sed -i '1i\RSID\tCHR\tPOS\tEffect_Allele\tOther_Allele\tBETA\tSE\tEffect_Allele_Frequency\tPVALUE\tN' \
sentinel_windows/liftover/age_at_menopause_hg38_with_samplesizes.txt

rm sentinel_windows/liftover/tmp_*
rm sentinel_windows/liftover/*.bed

# Get windows of sumstats (+/- 50kb) around each infertility lead SNP
for finf in female_infertility_analysis1 female_infertility_analysis3 female_infertility_analysis4 female_infertility_analysis5; do
	# Read lead SNPs
	while IFS=$'\t' read -r -a SNP_INFO; do
		awk -v CHR="${SNP_INFO[1]}" -v POS="${SNP_INFO[2]}" '$2=="chr"CHR && $3>=(POS-50000) && $3<=(POS+50000)' \
		sentinel_windows/liftover/age_at_menopause_hg38_with_samplesizes.txt \
		> sentinel_windows/${SNP_INFO[0]}_age_at_menopause.txt
		sed -i '1i\RSID\tCHR\tPOS\tEffect_Allele\tOther_Allele\tBETA\tSE\tEffect_Allele_Frequency\tPVALUE\tN' \
		sentinel_windows/${SNP_INFO[0]}_age_at_menopause.txt
	done < infertility_sentinels/${finf}_eur_sentinels.txt
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
