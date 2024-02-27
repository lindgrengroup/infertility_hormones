#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 2
#SBATCH -J filter_plot_UKBB_WES_results
#SBATCH -o /well/lindgren/samvida/hormones_infertility/exome_seq_results/logs/filter_plot_UKBB_WES_results-%j.out

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

cd /well/lindgren/samvida/hormones_infertility/exome_seq_results

# Merge gene-based results files across chromosomes
# head -n 1 results/1_${STRATA}_gene.txt > cat_results/${STRATA}_gene.txt
# tail -n +2 -q results/*_${STRATA}_gene.txt >> cat_results/${STRATA}_gene.txt

# Merge variant-based results files across chromosomes
# head -n 1 results/1_${STRATA}_variant.txt > cat_results/${STRATA}_variant.txt
# tail -n +2 -q results/*_${STRATA}_variant.txt >> cat_results/${STRATA}_variant.txt

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren/samvida/hormones_infertility/exome_seq_results/scripts/2_filter_plot_gene_results.R \
--strata=${STRATA}

Rscript /well/lindgren/samvida/hormones_infertility/exome_seq_results/scripts/2_filter_plot_variant_results.R \
--strata=${STRATA}


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
