#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J make_sumstats_gcta
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/make_sumstats_gcta-%j.out

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

# Load PLINK module
module load PLINK/2.00a2.3_x86_64

cd /well/lindgren/samvida/hormones_infertility/conditional_analysis/unreported_snps/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}

# Create temporary bed file of variants to extract
# chrN, pos0, pos1, varid
awk -F '\t' 'BEGIN {OFS=FS} NR>1 { print "chr"$2, $3-1, $3, $1}' \
/well/lindgren/samvida/hormones_infertility/meta_results_230613/filtered/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_filtered.txt \
> tmp_chrpos.bed

plink2 \
--bfile /well/lindgren/samvida/Resources/1000Genomes/hg38/all_hg38 \
--allow-extra-chr \
--extract range tmp_chrpos.bed \
--make-just-bim \
--threads 3 \
--memory 15000 \
--out tmp_chrpos_with1kg

# Merge the 1KG variant ID into summary statistics
awk '{ print "chr"$1":"$4, $2 }' \
tmp_chrpos_with1kg.bim | uniq | sort -k1 > tmp_chrpos_with1kg_variants.txt

join -1 1 -2 1 -o 1.1 2.2 1.5 1.6 1.7 1.9 1.10 1.11 \
<(sort -k1 /well/lindgren/samvida/hormones_infertility/meta_results_230613/filtered/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_filtered.txt) \
tmp_chrpos_with1kg_variants.txt \
> tmp_chrpos_gcta_sumstats.txt

awk -v N="$SAMPLE_SIZE" -F ' ' -v OFS="\t" '{ print $2,$3,$4,$5,$6,$7,$8,N}' \
tmp_chrpos_gcta_sumstats.txt > tmp_sumstats_gcta.txt
sed -i '1s/^/SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n/' tmp_sumstats_gcta.txt

rm tmp_chrpos*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
