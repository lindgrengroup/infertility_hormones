#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J filter_plot_public_results
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/filter_plot_public_results-%j.out
#SBATCH --array 1-17:1

echo "########################################################"
echo "Slurm Job ID: $SLURM_ARRAY_TASK_ID" 
echo "Slurm Task ID: $SLURM_ARRAY_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren/samvida/hormones_infertility/scripts/filter_plot_indiv_results.R \
--index=${SLURM_ARRAY_TASK_ID} 

# after liftOver
Rscript /well/lindgren/samvida/hormones_infertility/scripts/filter_plot_lifted_sumstats.R \
--index=${SLURM_ARRAY_TASK_ID} 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
