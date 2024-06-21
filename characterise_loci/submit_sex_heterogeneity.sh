#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J get_sex_heterogeneity_testosterone
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/get_sex_heterogeneity_testosterone-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

HORMONE=Testosterone

cd /well/lindgren/samvida/hormones_infertility

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript scripts/sex_heterogeneity.R \
--hormone=${HORMONE}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
