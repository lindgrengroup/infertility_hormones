#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J coloc_fi3_eur

#SBATCH --output /well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL/logs/coloc_fi3_eur-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL/logs/coloc_fi3_eur-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=4
#SBATCH --array=1-49



echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL

TRAIT="female_infertility_analysis3_eur"

config=/well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL/gtex_filelist.txt

# Extract the tissue name for the current $SLURM_ARRAY_TASK_ID
TISSUE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)



echo "This is array task ${SLURM_ARRAY_TASK_ID}, the tissue is ${TISSUE} and the trait is  ${TRAIT}." >> ${TRAIT}_output.txt

# Run coloc for TRAIT vs TISSUE
module load R/3.6.2-foss-2019b

Rscript run_coloc_infertility_eQTLs_gtex.R ${TRAIT} ${TISSUE}
