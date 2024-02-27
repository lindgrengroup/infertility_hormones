#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/01/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 10
#SBATCH -J regenie_step2_X
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/log_files/regenie_step2_X-%j.out

echo "########################################################"
echo "Slurm Job ID: ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

#HORMONE="FSH"
#SEX_STRATA="F"

mkdir /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step2

# if [ "${ANC_GROUP}" == "non_finnish_EUR" ]; then
# 	# X chromosome
# 	/well/lindgren/samvida/REGENIE/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
# 	--step 2 \
# 	--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chrX_v3.bgen \
# 	--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrX_v3_s486743.sample \
# 	--phenoFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/${HORMONE}_${SEX_STRATA}_cross_sectional.txt \
# 	--covarFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/sample_qc/qcd_sample_file_for_regenie.txt \
# 	--covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
# 	--catCovarList UKB_assmt_centre,sex,genotyping.array \
# 	--maxCatLevels 106 \
# 	--phenoColList ${HORMONE} \
# 	--qt \
# 	--pred /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step1/step1_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_pred.list \
# 	--bsize 1000 \
# 	--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step2/step2_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_chrX \
# 	--threads 10

# 	# XY chromosome (i.e. PAR region of Y)
# 	/well/lindgren/samvida/REGENIE/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
# 	--step 2 \
# 	--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chrXY_v3.bgen \
# 	--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrXY_v3_s486429.sample \
# 	--phenoFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/${HORMONE}_${SEX_STRATA}_cross_sectional.txt \
# 	--covarFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/sample_qc/qcd_sample_file_for_regenie.txt \
# 	--covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
# 	--catCovarList UKB_assmt_centre,sex,genotyping.array \
# 	--maxCatLevels 106 \
# 	--phenoColList ${HORMONE} \
# 	--qt \
# 	--pred /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step1/step1_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_pred.list \
# 	--bsize 1000 \
# 	--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step2/step2_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_chrXY \
# 	--threads 10
# else
# 	# X chromosome
# 	/well/lindgren/samvida/REGENIE/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
# 	--step 2 \
# 	--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chrX_v3.bgen \
# 	--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrX_v3_s486743.sample \
# 	--phenoFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_cross_sectional.txt \
# 	--covarFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/sample_qc/qcd_non_eur_sample_file_for_regenie.txt \
# 	--covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
# 	--catCovarList UKB_assmt_centre,sex,genotyping.array \
# 	--maxCatLevels 106 \
# 	--phenoColList ${HORMONE} \
# 	--qt \
# 	--pred /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step1/step1_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_pred.list \
# 	--bsize 1000 \
# 	--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step2/step2_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_chrX \
# 	--threads 10

# 	# XY chromosome (i.e. PAR region of Y)
# 	/well/lindgren/samvida/REGENIE/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
# 	--step 2 \
# 	--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chrXY_v3.bgen \
# 	--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrXY_v3_s486429.sample \
# 	--phenoFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_cross_sectional.txt \
# 	--covarFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/sample_qc/qcd_non_eur_sample_file_for_regenie.txt \
# 	--covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
# 	--catCovarList UKB_assmt_centre,sex,genotyping.array \
# 	--maxCatLevels 106 \
# 	--phenoColList ${HORMONE} \
# 	--qt \
# 	--pred /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step1/step1_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_pred.list \
# 	--bsize 1000 \
# 	--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step2/step2_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_chrXY \
# 	--threads 10
# fi

# X chromosome
/well/lindgren/samvida/REGENIE/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
--step 2 \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chrX_v3.bgen \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrX_v3_s486743.sample \
--phenoFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/${HORMONE}_${SEX_STRATA}_cross_sectional.txt \
--covarFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/sample_qc/qcd_sample_file_for_regenie.txt \
--covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
--catCovarList UKB_assmt_centre,sex,genotyping.array \
--maxCatLevels 106 \
--phenoColList ${HORMONE} \
--qt \
--pred /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step1/step1_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_pred.list \
--bsize 1000 \
--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step2/step2_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_chrX \
--threads 10

# XY chromosome (i.e. PAR region of Y)
/well/lindgren/samvida/REGENIE/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
--step 2 \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chrXY_v3.bgen \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrXY_v3_s486429.sample \
--phenoFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/traits_for_gwas/${HORMONE}_${SEX_STRATA}_cross_sectional.txt \
--covarFile /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/sample_qc/qcd_sample_file_for_regenie.txt \
--covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
--catCovarList UKB_assmt_centre,sex,genotyping.array \
--maxCatLevels 106 \
--phenoColList ${HORMONE} \
--qt \
--pred /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step1/step1_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_pred.list \
--bsize 1000 \
--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${HORMONE}_${ANC_GROUP}_${SEX_STRATA}/step2/step2_out_${HORMONE}_${ANC_GROUP}_${SEX_STRATA}_chrXY \
--threads 10


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
