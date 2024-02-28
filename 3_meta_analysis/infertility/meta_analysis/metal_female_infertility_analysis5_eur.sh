#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J metal_fi5_eur

#SBATCH --output /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/logs/metal_fi5_eur-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/logs/metal_fi5_eur-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=3



echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023

ANALYSIS="female_infertility_analysis5"
ANCESTRY="eur"

/apps/well/metal/20110325/metal metal_${ANALYSIS}_${ANCESTRY}.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

