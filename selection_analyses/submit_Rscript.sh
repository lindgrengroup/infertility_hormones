#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 15/03/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J get_gwascat_background_selection
#SBATCH -o /well/lindgren/samvida/hormones_infertility/selection_analyses/logs/get_gwascat_background_selection-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
Rscript /well/lindgren/samvida/hormones_infertility/scripts/tmp.R

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

