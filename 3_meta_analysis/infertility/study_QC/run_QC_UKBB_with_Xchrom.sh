#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J UKBB_QC

#SBATCH --output /well/lindgren/laura/projects/infertility/easyQC/logs/UKBB_QC-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/easyQC/logs/UKBB_QC-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=3



echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################
module load R/3.6.2-foss-2019b


for TRAIT in "female_infertility_analysis1" "female_infertility_analysis4" "male_infertility"
	do
	echo "1. Generate QQ plots"
	Rscript /well/lindgren/laura/projects/infertility/plots/qq_study_level/QQplots.R "UKBB" ${TRAIT}

	echo "2. Generate MH plots"
	Rscript /well/lindgren/laura/projects/infertility/plots/mh_study_level/make_MH_plots_results.R "UKBB" ${TRAIT}
	done 

for TRAIT in  "female_infertility_analysis2" "female_infertility_analysis3" "female_infertility_analysis5" 
	do
	echo "1. Generate QQ plots"
	Rscript /well/lindgren/laura/projects/infertility/plots/qq_study_level/QQplots_low_case_N.R "UKBB" ${TRAIT}
	echo "2. Generate MH plots"
	Rscript /well/lindgren/laura/projects/infertility/plots/mh_study_level/make_MH_plots_results.R "UKBB" ${TRAIT}
	done 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

