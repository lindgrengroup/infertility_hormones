#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J grep_exposures_outcomes
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/twosmr_grep_exposures_outcomes-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/two_sample_mr

# Add RSIDs to infertility files
# Convert chr:pos to RSIDs
while IFS=$'\t' read -r -a STRATA_INFO
do
	# Create a chr:pos column, make sure the position is an integer, 
	# make the alleles capital letters, make the sample sizes integers,
	# subset to only the columns that are needed for MR, and sort by the chr:pos column
	awk -v OFS="\t" '{print "chr"$12":"($13+0), $12, ($13+0), toupper($2), toupper($3), $4, $5, $6, $7, ($9+0), ($10+0)}' \
	${STRATA_INFO[1]} | sort -k1 > outcomes/tmp_${STRATA_INFO[0]}_sumstats.txt

	# Join by chr:pos
	join -t $'\t' -1 1 -2 6 -o 2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.3,2.4 \
	outcomes/tmp_${STRATA_INFO[0]}_sumstats.txt \
	/well/lindgren/samvida/Resources/hg38/grch38.withchrpos.all.rsid.txt \
	> outcomes/tmp_${STRATA_INFO[0]}_sumstats_with_rsids.txt

	# Remove rows where the map alleles don't match the sumstat alleles
	awk 'BEGIN { FS = OFS = "\t" } { split($13, values, ","); for (i in values) if ($4 == $12 && $5 == values[i] || $5 == $12 && $4 == values[i]) { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11; break } }' \
	outcomes/tmp_${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	> outcomes/${STRATA_INFO[0]}_sumstats_with_rsids.txt

	# Add header column
	sed -i '1s/^/RSID\tCHR\tPOS\tALLELE1\tALLELE2\tFREQ1\tBETA\tSE\tPVALUE\tN_CASES\tN_CONTROLS\n/' \
	outcomes/${STRATA_INFO[0]}_sumstats_with_rsids.txt

	rm outcomes/tmp_*

done < infertility_EUR_files_list.txt

# Concatenate all the exposure RSIDs into one file
awk 'NR > 1 {print $2}' exposures/all_lead_snp_sumstats_*_with_rsids.txt | sort | uniq > exposures/all_hormones_rsids.txt

# Look up these variants in infertility results
while IFS=$'\t' read -r -a STRATA_INFO
do
	grep -wf exposures/all_hormones_rsids.txt outcomes/${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	> outcomes/instrument_outcomes_for_${STRATA_INFO[0]}.txt
done < infertility_EUR_files_list.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

