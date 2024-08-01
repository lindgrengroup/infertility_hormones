#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J metal_fi1_all

#SBATCH --output /well/lindgren/samvida/hormones_infertility/infertility_meta_mvp/logs/metal_fi1_all-%j.out 

#SBATCH --error /well/lindgren/samvida/hormones_infertility/infertility_meta_mvp/logs/metal_fi1_all-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=3

echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/infertility_meta_mvp

ANALYSIS="female_infertility_analysis1"
ANCESTRY="all"

/apps/well/metal/20110325/metal scripts/metal_${ANALYSIS}_${ANCESTRY}.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

