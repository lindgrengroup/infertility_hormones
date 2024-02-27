#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J munge_sumstats_ldsc
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/munge_sumstats_ldsc-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility
mkdir munged_sumstats

# Convert chr:pos to RSIDs
# To speed this up, create a version of the hg38 chrpos map file that is limited to HapMap3 IDs since
# we will need to subset to these later anyway

# awk 'FNR==NR { hm3_rsids[$1]; next } $5 in hm3_rsids' \
# /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
# /well/lindgren/samvida/Resources/hg38/grch38.withchrpos.all.rsid.txt \
# | sort -k6 > /well/lindgren/samvida/Resources/hg38/grch38.withchrpos.hm3.rsid.txt

# For hormones
while IFS=$'\t' read -r -a STRATA_INFO
do
	# Join by chr:pos
	join -t $'\t' -1 1 -2 6 -o 2.5,1.5,1.6,1.7,1.9,1.11,2.3,2.4 \
	<(sort -k1 ${STRATA_INFO[1]}) \
	/well/lindgren/samvida/Resources/hg38/grch38.withchrpos.hm3.rsid.txt \
	> munged_sumstats/tmp_${STRATA_INFO[0]}_sumstats_with_rsids.txt

	# Remove rows where the map alleles don't match the sumstat alleles
	awk 'BEGIN { FS = OFS = "\t" } { split($8, values, ","); for (i in values) if ($2 == $7 && $3 == values[i] || $3 == $7 && $2 == values[i]) { print $1,$2,$3,$4,$5,$6; break } }' \
	munged_sumstats/tmp_${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	> munged_sumstats/${STRATA_INFO[0]}_sumstats_with_rsids.txt

	# Add header column
	sed -i '1s/^/RSID\tAllele1\tAllele2\tFreq1\tBETA\tPVALUE\n/' \
	munged_sumstats/${STRATA_INFO[0]}_sumstats_with_rsids.txt

	rm munged_sumstats/tmp_*

done < hormone_files_list.txt

# For infertility traits
while IFS=$'\t' read -r -a STRATA_INFO
do
	# Create a chr:pos column, make sure the position is an integer, 
	# make the alleles capital letters, make the sample sizes integers,
	# subset to only the columns that are needed for LDSC, and sort by the chr:pos column
	awk -v OFS="\t" '{print "chr"$12":"($13+0), toupper($2), toupper($3), $4, $5, $7, ($9+0), ($10+0)}' \
	${STRATA_INFO[1]} | sort -k1 > munged_sumstats/tmp_${STRATA_INFO[0]}_sumstats.txt

	# Join by chr:pos
	join -t $'\t' -1 1 -2 6 -o 2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.3,2.4 \
	munged_sumstats/tmp_${STRATA_INFO[0]}_sumstats.txt \
	/well/lindgren/samvida/Resources/hg38/grch38.withchrpos.hm3.rsid.txt \
	> munged_sumstats/tmp_${STRATA_INFO[0]}_sumstats_with_rsids.txt

	# Remove rows where the map alleles don't match the sumstat alleles
	awk 'BEGIN { FS = OFS = "\t" } { split($10, values, ","); for (i in values) if ($2 == $9 && $3 == values[i] || $3 == $9 && $2 == values[i]) { print $1,$2,$3,$4,$5,$6,$7,$8; break } }' \
	munged_sumstats/tmp_${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	> munged_sumstats/${STRATA_INFO[0]}_sumstats_with_rsids.txt

	# Add header column
	sed -i '1s/^/RSID\tAllele1\tAllele2\tFreq1\tEffect\tP.value\tN_CASES\tN_CONTROLS\n/' \
	munged_sumstats/${STRATA_INFO[0]}_sumstats_with_rsids.txt

	rm munged_sumstats/tmp_*

done < infertility_files_list.txt

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

# For hormone traits
while IFS=$'\t' read -r -a STRATA_INFO
do
	munge_sumstats.py \
	--sumstats munged_sumstats/${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	--N ${STRATA_INFO[2]} \
	--snp RSID --a1 Allele1 --a2 Allele2 --p PVALUE \
	--frq Freq1 --signed-sumstats BETA,0 \
	--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
	--out munged_sumstats/${STRATA_INFO[0]}_for_ldsc
done < hormone_files_list.txt

# For infertility traits
while IFS=$'\t' read -r -a STRATA_INFO
do
	munge_sumstats.py \
	--sumstats munged_sumstats/${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	--N-cas-col N_CASES \
	--N-con-col N_CONTROLS \
	--snp RSID --a1 Allele1 --a2 Allele2 --p P.value \
	--frq Freq1 --signed-sumstats Effect,0 \
	--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
	--out munged_sumstats/${STRATA_INFO[0]}_for_ldsc
done < infertility_files_list.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0



