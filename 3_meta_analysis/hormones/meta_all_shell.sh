#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/01/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J meta_analysis_no_ukb_230613
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/meta_analysis_no_ukb_230613-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

MAIN_FILEPATH="/well/lindgren/samvida/hormones_infertility"
cd $MAIN_FILEPATH

# cat strata_list.txt | while read STRATA;
# do
# 	if test -f ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/meta_${STRATA}_EUR_script.txt; then
# 		/apps/well/metal/20110325/metal < ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/meta_${STRATA}_EUR_script.txt
# 	fi

# 	if test -f ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/meta_${STRATA}_all_anc_script.txt; then
# 		/apps/well/metal/20110325/metal < ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/meta_${STRATA}_all_anc_script.txt
# 	fi
# done

cat strata_list.txt | while read STRATA;
do
	if test -f ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/no_UKB/meta_${STRATA}_EUR_script.txt; then
		/apps/well/metal/20110325/metal < ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/no_UKB/meta_${STRATA}_EUR_script.txt
	fi

	if test -f ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/no_UKB/meta_${STRATA}_all_anc_script.txt; then
		/apps/well/metal/20110325/metal < ${MAIN_FILEPATH}/scripts/meta_analysis_instructions/no_UKB/meta_${STRATA}_all_anc_script.txt
	fi
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
