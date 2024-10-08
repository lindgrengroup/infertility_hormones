#!/bin/bash

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH --array 2,4,6,7,11,12,16,17,19,20
#SBATCH -J extract_genotypes

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load PLINK/2.00a2.3_x86_64

CHR_GET=${SLURM_ARRAY_TASK_ID}

cd /well/lindgren/samvida/hormones_infertility/exome_seq_results/conditioned_on_common

awk -v chr="${CHR_GET}" '$2 == chr {print $1}' nearby_common_variants/unique_rsids_to_get.txt \
> tmp_chr${CHR_GET}_rsids.txt

plink2 \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${CHR_GET}_v3.bgen \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--extract tmp_chr${CHR_GET}_rsids.txt \
--make-pgen \
--out plink_files/chr${CHR_GET}_rsids \
--memory 15000 \
--threads 5 

rm tmp_chr${CHR_GET}_rsids.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
