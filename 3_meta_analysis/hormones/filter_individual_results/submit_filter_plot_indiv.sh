#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 2
#SBATCH -J filter_plot_ALSPAC_results
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/filter_plot_ALSPAC_results-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren/samvida/hormones_infertility/scripts/filter_plot_indiv_results_BOLT.R \
--inputFile=${inputFile} \
--logFile=${logFile} \
--outputFile=${outputFile} \
--outPlotDir=${outPlotDir} 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
