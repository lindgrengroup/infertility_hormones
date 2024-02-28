#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J step2_regenie_chrX

#SBATCH --output /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/scripts/eur/logs/step2_regenie_chrX-%j.out 

#SBATCH --error /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/scripts/eur/logs/step2_regenie_chrX-%j.err 

#SBATCH -p long

#SBATCH --cpus-per-task=10





echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load BGEN/1.1.6-GCCcore-8.3.0
​

/well/lindgren/UKBIOBANK/laura/REGENIE/regenie_v3.2.2.gz_x86_64_Centos7_mkl \
      --step 2 \
      --bgen /well/ukbb-wtchg/v3/imputation/ukb_imp_chrX_v3.bgen \
      --sample /well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chrX_v3_s486743.sample \
      --ref-first \
      --phenoFile /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/eur_only/regenie_pheno_input_file_all_traits_EUR_ONLY.txt \
	--phenoColList female_infertility,anovulatory_infertility,anatomical_infertility,idiop_infertility_inclusion,idiop_infertility_exclusion,male_infertility \
      --covarFile /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/eur_only/regenie_covar_input_file_all_traits_EUR_ONLY.txt \
      --covarColList assessment_centre,sex,genotyping.array,batch,age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
      --catCovarList assessment_centre,sex,genotyping.array,batch \
      --maxCatLevels 106 \
      --bt \
      --firth --approx --firth-se \
      --pred /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/output/europeans/step1_out_pred.list \
      --bsize 400 \
      --niter 100 \
      --gz \
      --threads 22 \
      --minINFO 0.8 \
      --out /well/lindgren/UKBIOBANK/laura/infertility/REGENIE/output/europeans/step2_out_chrX \

​
echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0















