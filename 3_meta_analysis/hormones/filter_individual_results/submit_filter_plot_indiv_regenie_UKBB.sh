#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 2
#SBATCH -J filter_plot_UKBB_results
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/filter_plot_UKBB_results-%j.out

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

# Cat results and move the summary statistics to the full data folder
cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/GWAS/REGENIE_results/${PHENO}_${ANC_GROUP}_${SS}/step2
# Concatenate results across chromosomes
head -n 1 step2_out_${PHENO}_${ANC_GROUP}_${SS}_chr1_${PHENO}.regenie > ${PHENO}_${ANC_GROUP}_${SS}.regenie
tail -n +2 -q step2_out_${PHENO}_${ANC_GROUP}_${SS}_chr*_${PHENO}.regenie >> ${PHENO}_${ANC_GROUP}_${SS}.regenie

cp ${PHENO}_${ANC_GROUP}_${SS}.regenie /well/lindgren/samvida/hormones_infertility/data/UKBB

MAIN_FILEPATH="/well/lindgren/samvida/hormones_infertility"
cd ${MAIN_FILEPATH}/data/UKBB

# Create log file for filtering results
LOG_FILE="${MAIN_FILEPATH}/data/UKBB/filtered/${PHENO}_${ANC_GROUP}_${SS}.log"
rm $LOG_FILE
touch $LOG_FILE

# Create temporary directory for filtering
rm -r ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC
mkdir ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC

PREQC=$(wc -l < ${PHENO}_${ANC_GROUP}_${SS}.regenie)
printf "** Phenotype, ancestry, strata: ${PHENO}_${ANC_GROUP}_${SS}\n" >> $LOG_FILE 
printf "\t # SNPs pre-QC: $((${PREQC}-1)) \n" >> $LOG_FILE

# Location of filtering files
FILTER_LOC="/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/ukb_imp_snp_qc"

# INFO filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $3 in a' \
${FILTER_LOC}/passed_info_QC.txt ${PHENO}_${ANC_GROUP}_${SS}.regenie \
> ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_info.txt

TMPQC=$(wc -l < ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_info.txt)
printf "\t # SNPs passed INFO > 0.8: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Genotyping filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $3 in a' \
${FILTER_LOC}/${ANC_GROUP}/passed_geno_QC.txt \
./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_info.txt \
> ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_geno.txt

TMPQC=$(wc -l < ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_geno.txt)
printf "\t # SNPs passed missingness < 0.05: $((${TMPQC}-1)) \n" >> $LOG_FILE

# HWE filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $3 in a' \
${FILTER_LOC}/${ANC_GROUP}/passed_hwe_QC.txt \
./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_geno.txt \
> ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_hwe.txt

TMPQC=$(wc -l < ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_hwe.txt)
printf "\t # SNPs passed HWE pval > 1E-06: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Biallelic filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $3 in a' \
${FILTER_LOC}/${ANC_GROUP}/passed_biallelic_QC.txt \
./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_hwe.txt \
> ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_biallelic.txt

TMPQC=$(wc -l < ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_biallelic.txt)
printf "\t # SNPs passed bi-allelic QC: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Save results
cp ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC/tmp_passed_biallelic.txt ${PHENO}_${ANC_GROUP}_${SS}_init_qc.txt
rm -r ./${PHENO}_${ANC_GROUP}_${SS}_tmp_QC

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren/samvida/hormones_infertility/scripts/filter_plot_indiv_results_regenie.R \
--inputFile=${MAIN_FILEPATH}/data/UKBB/${PHENO}_${ANC_GROUP}_${SS}_init_qc.txt \
--logFile=$LOG_FILE \
--outputFile=${MAIN_FILEPATH}/data/UKBB/filtered/${PHENO}_${ANC_GROUP}_${SS}_filtered.txt \
--outPlotDir=${MAIN_FILEPATH}/plots/UKBB/${PHENO}_${ANC_GROUP}_${SS}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
