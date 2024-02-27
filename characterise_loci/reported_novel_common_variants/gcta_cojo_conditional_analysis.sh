#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J gcta_cojo_conditional_analysis
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/gcta_cojo_conditional_analysis-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

# Load PLINK module
module load PLINK/2.00a2.3_x86_64

cd /well/lindgren/samvida/hormones_infertility/conditional_analysis/unreported_snps/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}

mkdir tmp_1000G_bed
mkdir tmp_1000G_ld

touch ${HORMONE_TYPE}_reported_snp_list.txt
touch ${HORMONE_TYPE}_poss_novel_snp_list.txt

while read var; do
	# Create a temporary file for single variant	
	echo ${var} > tmp_${var}.txt
	# Remove SNP of interest from published SNP list
	grep -vwE ${var} \
	${HORMONE_TYPE}_published_snps_around_${var}.txt > tmp_${HORMONE_TYPE}_published_snps_around_${var}.txt
			
	# Create temporary bed/bim/fam files for the relevant genomic region
	# Keep samples of just EUR ancestry or all ancestries based on the group
	if [[ "${ANC_GROUP}" == "EUR" ]]; then
		plink2 \
		--bfile /well/lindgren/samvida/Resources/1000Genomes/hg38/all_hg38 \
		--allow-extra-chr \
		--keep /well/lindgren/samvida/Resources/1000Genomes/hg38/sample_iids_EUR.txt \
		--snp ${var} \
		--window 1000 \
		--rm-dup force-first \
		--make-bed \
		--threads 3 \
		--memory 15000 \
		--out tmp_1000G_bed/${var}_region
	else
		plink2 \
		--bfile /well/lindgren/samvida/Resources/1000Genomes/hg38/all_hg38 \
		--allow-extra-chr \
		--snp ${var} \
		--window 1000 \
		--rm-dup force-first \
		--make-bed \
		--threads 3 \
		--memory 15000 \
		--out tmp_1000G_bed/${var}_region
	fi

	# Create temporary LD-matrix for SNP of interest vs published SNPs in the region
	/apps/well/plink/1.90b3/plink \
	--bfile tmp_1000G_bed/${var}_region \
	--r2 inter-chr \
	--ld-snp ${var} \
	--threads 3 \
	--memory 15000 \
	--out tmp_1000G_ld/tmp_${var}_region
			
	# Only retain SNPs of interest (published)
	grep -f tmp_${HORMONE_TYPE}_published_snps_around_${var}.txt \
	tmp_1000G_ld/tmp_${var}_region.ld > tmp_1000G_ld/${var}_region.ld
	rm tmp_1000G_ld/tmp_${var}_region*

	### CONDITIONAL ANALYSIS

	# First remove collinear SNPs from conditional SNP-list
	/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 \
	--bfile tmp_1000G_bed/${var}_region  \
	--cojo-file tmp_sumstats_gcta.txt \
	--extract ${HORMONE_TYPE}_published_snps_around_${var}.txt \
	--cojo-slct \
	--out ${var}_indpt_published_snps_${HORMONE_TYPE}

	# If there are no independent SNPs then the SNP of interest
	# is previously reported
	if [[ ! -e "${var}_indpt_published_snps_${HORMONE_TYPE}.jma.cojo" ]]; then
		echo ${var} >> ${HORMONE_TYPE}_reported_snp_list.txt
	else
		# Extract independent SNP list
		awk -F '\t' '(NR>1) {print $2}' \
		${var}_indpt_published_snps_${HORMONE_TYPE}.jma.cojo \
		> ${var}_indpt_published_snps_${HORMONE_TYPE}.txt

		# If our SNP of interest is in the independently associated 
		# list of SNPs, then it may be novel;
		# else it is reported
		if grep -q ${var} ${var}_indpt_published_snps_${HORMONE_TYPE}.txt; then 
			# Remove SNP of interest from reported SNP list
			grep -vwE ${var} ${var}_indpt_published_snps_${HORMONE_TYPE}.txt \
			> ${var}_indpt_published_snps_${HORMONE_TYPE}_for_cond.txt	

			# If there are other SNPs in the independently associated list,
			# our SNP may or may not be novel, have to check conditional associations
			# If this was the only SNP that remained, it is perfectly collinear so it is reported
			if [ -s ${var}_indpt_published_snps_${HORMONE_TYPE}_for_cond.txt ]; then
				echo ${var} >> ${HORMONE_TYPE}_poss_novel_snp_list.txt
				# Get effect of SNP of interest, conditioned on remaining SNPs
				/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 \
				--bfile tmp_1000G_bed/${var}_region  \
				--cojo-file tmp_sumstats_gcta.txt \
				--extract ${var}_indpt_published_snps_${HORMONE_TYPE}.txt \
				--cojo-cond ${var}_indpt_published_snps_${HORMONE_TYPE}_for_cond.txt \
				--out ${var}_cond_on_pub_${HORMONE_TYPE}
			else
				echo ${var} >> ${HORMONE_TYPE}_reported_snp_list.txt
			fi
		else
			echo ${var} >> ${HORMONE_TYPE}_reported_snp_list.txt
		fi
	fi
done < unreported_${HORMONE_TYPE}.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

