#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J add_rsids_hormones
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/add_rsids_hormones-%j.out

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

cd /well/lindgren/samvida/hormones_infertility/meta_results_230613

mkdir add_rsids

# Sort by chr:position to add rsid
# Replace "chr23" with "chrX" in order to get the rsids
sed 's/chr23/chrX/g' filtered/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_filtered.txt \
> add_rsids/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_filtered.txt

join -t $'\t' -1 1 -2 6 -o 1.1,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.3,2.4 \
<(sort -k1 add_rsids/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_filtered.txt) \
/well/lindgren/samvida/Resources/hg38/grch38.withchrpos.all.rsid.txt \
> add_rsids/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_tmp.txt

# Remove rows where the map alleles don't match the sumstat alleles
awk 'BEGIN { FS = OFS = "\t" } { split($16, values, ","); for (i in values) if ($6 == $15 && $7 == values[i] || $7 == $15 && $6 == values[i]) { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; break } }' \
add_rsids/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_tmp.txt \
> add_rsids/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_with_rsids.txt

# Add header column
sed -i '1s/^/ID\tRSID\tCHROM\tGENPOS\tMAF\tAllele1\tAllele2\tFreq1\tFreqSE\tBETA\tSE\tPVALUE\tDirection\tHetPVal\n/' \
add_rsids/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_with_rsids.txt

rm add_rsids/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_tmp.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
