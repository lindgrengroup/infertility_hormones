#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J prep_coloc_sumstats_inhbb
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/prep_coloc_sumstats_inhbb-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility

# Lift AMH sumstats up to hg38 and add sample sizes
module load liftOver/20210519

# Create a chr:pos column, make the alleles letters removing quotes, 
# and sort by the chr:pos column
dos2unix rg_repro_traits/external_sumstats/verdiesen_AMH_231127_buildGRCh37.tsv # to get rid of the file-endings
awk -v OFS="\t" '{gsub(/"/, "", $4); gsub(/"/, "", $5); print "chr"$1":"$2, $1, $2, $4, $5, $6, $7, $8, $3}' \
rg_repro_traits/external_sumstats/verdiesen_AMH_231127_buildGRCh37.tsv \
| sort -k1 > colocalisation/sentinel_windows/liftover/tmp_verdiesen_AMH_sumstats.txt

join -t $'\t' -1 1 -2 6 -o 2.3,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.4,2.5 \
colocalisation/sentinel_windows/liftover/tmp_verdiesen_AMH_sumstats.txt \
<(sort -k6 /well/lindgren/samvida/Resources/hg19/grch37.withchrpos.hm3.rsid.txt) \
> colocalisation/sentinel_windows/liftover/tmp_verdiesen_AMH_sumstats_with_rsids.txt

# Remove rows where the map alleles don't match the sumstat alleles
awk 'BEGIN { OFS = "\t" } { split($11, values, ","); for (i in values) if ($4 == $10 && $5 == values[i] || $5 == $10 && $4 == values[i]) { print $1,$2,$3,$4,$5,$6,$7,$8,$9; break } }' \
colocalisation/sentinel_windows/liftover/tmp_verdiesen_AMH_sumstats_with_rsids.txt \
> colocalisation/sentinel_windows/liftover/verdiesen_AMH_sumstats_with_rsids.txt

# Add header column
sed -i '1s/^/RSID\tCHR\tPOS\tAllele1\tAllele2\tFreq1\tBETA\tSE\tPVALUE\n/' \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_sumstats_with_rsids.txt
rm colocalisation/sentinel_windows/liftover/tmp_*.txt
	
# Create bed file with chrX, pos0, pos1, RSID
awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$2, $3-1, $3, $1}' \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_sumstats_with_rsids.txt > \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg19.bed
			
liftOver colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg19.bed \
/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg38.bed \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_unlifted.bed

# Create a new sumstats file by joining via the RSID and add sample sizes for N
# unique_id: 4th column in the lifted over bed file and 1st column in the original sumstats
join -1 4 -2 1 -o 2.1,1.1,1.3,2.4,2.5,2.6,2.7,2.8,2.9 \
<(sort -k 4 colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg38.bed) \
<(sort -k 1 colocalisation/sentinel_windows/liftover/verdiesen_AMH_sumstats_with_rsids.txt) \
> colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg38.txt
awk -v N="7049" 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, N }' \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg38.txt \
> colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg38_with_samplesizes.txt

sed -i '1i\RSID\tCHR\tPOS\tEffect_Allele\tOther_Allele\tEffect_Allele_Frequency\tBETA\tSE\tPVALUE\tN' \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg38_with_samplesizes.txt

rm colocalisation/sentinel_windows/liftover/tmp_*
rm colocalisation/sentinel_windows/liftover/*.bed

# Get windows of sumstats (+/- 50kb) around the INHBB SNP of interest
awk -v CHR="2" -v POS="120388925" '$2=="chr"CHR && $3>=(POS-50000) && $3<=(POS+50000)' \
colocalisation/sentinel_windows/liftover/verdiesen_AMH_hg38_with_samplesizes.txt \
> colocalisation/sentinel_windows/chr2_120388925_verdiesen_AMH.txt
sed -i '1i\RSID\tCHR\tPOS\tEffect_Allele\tOther_Allele\tBETA\tSE\tEffect_Allele_Frequency\tPVALUE\tN' \
colocalisation/sentinel_windows/chr2_120388925_verdiesen_AMH.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
