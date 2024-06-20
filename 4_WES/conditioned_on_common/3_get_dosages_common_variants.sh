#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 20/06/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J get_variant_dosages
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/get_variant_dosages-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/exome_seq_results/conditioned_on_common/dosage_files

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

UKB_PATH="/well/lindgren-ukbb/projects/ukbb-11867/DATA"

/apps/well/qctool/2.0.1/qctool \
-g ${UKB_PATH}/IMPUTATION/ukb_imp_chr${CHR}_v3.bgen \
-s ${UKB_PATH}/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
-condition-on ${VARID} \
-os ${VARID}_dosages.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
