#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/01/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH --array 1-22%4 
#SBATCH -J lava_bivariate
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/lava_bivariate-%j.out

echo "########################################################"
echo "Slurm Job ID: ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren/samvida/hormones_infertility/scripts/2_run_lava_bivariate.R \
--pheno1=${pheno1} \
--pheno2=${pheno2} \
--chrom=${SLURM_ARRAY_TASK_ID}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
