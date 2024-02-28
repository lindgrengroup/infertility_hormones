#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J sentinel

#SBATCH --output /well/lindgren/laura/projects/infertility/sentinel_snps/logs/sentinel-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/sentinel_snps/logs/sentinel-%j.err 
#SBATCH -p short

#SBATCH --cpus-per-task=3



echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################
module load R/3.6.2-foss-2019b


for TRAIT in "female_infertility_analysis1_all" "female_infertility_analysis1_eur" "female_infertility_analysis2_all" "female_infertility_analysis2_eur" "female_infertility_analysis3_all" "female_infertility_analysis3_eur" "female_infertility_analysis4_all" "female_infertility_analysis4_eur" "female_infertility_analysis5_all" "female_infertility_analysis5_eur" "male_infertility_all" "male_infertility_eur" ;
	do 
	Rscript /well/lindgren/laura/projects/infertility/sentinel_snps/identify_sentinel_snps.R ${TRAIT}
	done 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

