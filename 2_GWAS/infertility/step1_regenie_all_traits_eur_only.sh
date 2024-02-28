#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J step1_regenie

#SBATCH --output /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/scripts/eur/logs/step1_regenie-%j.out 

#SBATCH --error /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/scripts/eur/logs/step1_regenie-%j.err 

#SBATCH -p short 
#SBATCH --constraint=skl-compat 

#SBATCH --cpus-per-task=10




echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

/well/lindgren/UKBIOBANK/laura/REGENIE/regenie_v3.2.2.gz_x86_64_Centos7_mkl --step 1 \
--bed /well/lindgren/UKBIOBANK/laura/REGENIE/ld_clump/ukb_cal_v2_qced_pruned \
--phenoFile /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/eur_only/regenie_pheno_input_file_all_traits_EUR_ONLY.txt \
--covarFile /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/eur_only/regenie_covar_input_file_all_traits_EUR_ONLY.txt \
--covarColList assessment_centre,sex,genotyping.array,batch,age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
--catCovarList assessment_centre,sex,genotyping.array,batch \
--maxCatLevels 106 \
--phenoColList female_infertility,anovulatory_infertility,anatomical_infertility,idiop_infertility_inclusion,idiop_infertility_exclusion,IHH,IHH_female,IHH_male,POI,POI_by_aam,male_infertility \
--bt \
--bsize 1000 --out /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/output/europeans/step1_out \
--loocv  \
--threads 10

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

