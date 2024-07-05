#!/bin/bash

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J extract_genotypes_chrX

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load PLINK/2.00a2.3_x86_64

cd /well/lindgren/samvida/hormones_infertility/exome_seq_results/conditioned_on_common

awk '$2 == 23 {print $1}' nearby_common_variants/unique_rsids_to_get.txt \
> tmp_chrX_rsids.txt

plink2 \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chrX_v3.bgen \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrX_v3_s486743.sample \
--extract tmp_chrX_rsids.txt \
--make-pgen \
--out plink_files/chrX_rsids \
--memory 15000 \
--threads 5 

rm tmp_chrX_rsids.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
