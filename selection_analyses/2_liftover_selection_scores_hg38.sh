#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 02/10/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J get_selection_sumstats_infertility_loci
#SBATCH -o /well/lindgren/samvida/hormones_infertility/selection_analyses/logs/get_selection_sumstats_infertility_loci-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/Resources/selection_scores_public

# Liftover Beta2 and SDS scores from hg19 to hg38 positions
module load liftOver/20210519

# SDS scores
mkdir SDS_UK10K_n3195_release_Sep_19_2016/hg38

# Create bed file with chrX, pos0, pos1, unique_id (RSID)
awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$1, $2-1, $2+0, $3}' \
SDS_UK10K_n3195_release_Sep_19_2016/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals.txt \
> SDS_UK10K_n3195_release_Sep_19_2016/SDS_UK10K_n3195_release_Sep_19_2016_hg19.bed
			
liftOver SDS_UK10K_n3195_release_Sep_19_2016/SDS_UK10K_n3195_release_Sep_19_2016_hg19.bed \
/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_hg38.bed \
SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_unlifted.bed

# Merge new chromosome position with old file to get scores
# Create a new sumstats file by joining via the RSID: 4th column in the lifted over bed file and 3rd column in the original sumstats
# Add header
join -1 4 -2 3 -o 1.1,1.3,2.3,2.4,2.5,2.6,2.7,2.8 \
<(sort -k 4 SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_hg38.bed) \
<(sort -k 3 SDS_UK10K_n3195_release_Sep_19_2016/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals.txt) \
> SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals_hg38.txt
sed -i '1s/^/CHR POS ID AA DA DAF SDS pval\n/' \
SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals_hg38.txt

# Beta2 scores
mkdir BetaScan2/hg38
for CHR in {1..22}; do
	# Add unique chr_pos ID 
	awk -v FS="\t" -v OFS="\t" -v CHR="$CHR" '{print $0, "chr"CHR"_"$1}' \
	BetaScan2/StdB2Scores/chr${CHR}_B2std.out \
	> BetaScan2/StdB2Scores/chr${CHR}_B2std_with_unique_id.txt
	
	# Create bed file with chrX, pos0, pos1, unique_id (CHR_POS)
	awk -v FS="\t" -v OFS="\t" -v CHR="$CHR" 'NR > 1 {print "chr"CHR, $1-1, $1+0, $4}' \
	BetaScan2/StdB2Scores/chr${CHR}_B2std_with_unique_id.txt \
	> BetaScan2/StdB2Scores/chr${CHR}_B2std_hg19.bed
				
	liftOver BetaScan2/StdB2Scores/chr${CHR}_B2std_hg19.bed \
	/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
	BetaScan2/hg38/chr${CHR}_B2std_hg38.bed \
	BetaScan2/hg38/chr${CHR}_B2std_unlifted.bed

	# Merge new chromosome position with old file to get scores
	# Create a new sumstats file by joining via the RSID: 4th column in the lifted over bed file and 4th column in the original sumstats
	join -1 4 -2 4 -o 1.1,1.3,2.2,2.3 \
	<(sort -k 4 BetaScan2/hg38/chr${CHR}_B2std_hg38.bed) \
	<(sort -k 4 BetaScan2/StdB2Scores/chr${CHR}_B2std_with_unique_id.txt) \
	> BetaScan2/hg38/chr${CHR}_B2std_hg38.txt
done 
# Merge all chromosomes
cat BetaScan2/hg38/chr*_B2std_hg38.txt >> BetaScan2/hg38/all_B2std_hg38.txt
# Add header
sed -i '1s/^/Chromosome Position Beta2 Beta2_std\n/' \
BetaScan2/hg38/all_B2std_hg38.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0