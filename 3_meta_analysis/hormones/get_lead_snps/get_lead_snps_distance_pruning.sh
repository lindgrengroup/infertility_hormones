#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J get_lead_snps_distance_pruning
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/get_lead_snps_distance_pruning-%j.out

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

mkdir lead_snps

awk -F '\t' 'BEGIN {OFS = FS} {if ($11 <= 5E-08) print}' \
filtered/${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_filtered.txt \
> lead_snps/gws_variants_${HORMONE}_${SEX_STRATA}_${ANC_GROUP}.txt

cd lead_snps
# Sort by chr:position to add rsid
# Replace "chr23" with "chrX" in order to get the rsids
sed -i 's/chr23/chrX/g' gws_variants_${HORMONE}_${SEX_STRATA}_${ANC_GROUP}.txt

join -t $'\t' -1 1 -2 6 -o 1.1,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.3,2.4 \
<(sort -k1 gws_variants_${HORMONE}_${SEX_STRATA}_${ANC_GROUP}.txt) \
/well/lindgren/samvida/Resources/hg38/grch38.withchrpos.all.rsid.txt \
> ${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_tmp.txt

# Remove rows where the map alleles don't match the sumstat alleles
awk 'BEGIN { FS = OFS = "\t" } { split($16, values, ","); for (i in values) if ($6 == $15 && $7 == values[i] || $7 == $15 && $6 == values[i]) { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; break } }' \
${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_tmp.txt \
> gws_variants_${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_with_rsids.txt

# Add header column
sed -i '1s/^/ID\tRSID\tCHROM\tGENPOS\tMAF\tAllele1\tAllele2\tFreq1\tFreqSE\tBETA\tSE\tPVALUE\tDirection\tHetPVal\n/' \
gws_variants_${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_with_rsids.txt

rm gws_variants_${HORMONE}_${SEX_STRATA}_${ANC_GROUP}.txt
rm ${HORMONE}_${SEX_STRATA}_${ANC_GROUP}_tmp.txt

# Also create a subset of GWS variants that are present in the 1000 Genomes project (higher confidence + for classification)
# Load PLINK module
module load PLINK/2.00a2.3_x86_64

mkdir plink_log
# Read all GWS SNPs and create list of variants to check for
# Remvoe any SNPs on X chromosome because we can't classify these with LD anyway
awk '$3 != 23 {print $3":"$4":"$6":"$7}' gws_variants_*_with_rsids.txt | sort | uniq > plink_log/all_gws_varids.txt

# Extract bim file for variants of interest 
plink2 \
--bfile /well/lindgren/samvida/Resources/1000Genomes/hg38/all_hg38 \
--allow-extra-chr \
--extract plink_log/all_gws_varids.txt \
--make-just-bim \
--threads 3 \
--memory 15000 \
--out plink_log/checking_gws_varids

# Intersect the variants in the bim file with the original list to figure out which ones did not make it
awk '{print $2}' plink_log/checking_gws_varids.bim > plink_log/do_not_flip.txt
grep -v -F -x -f plink_log/do_not_flip.txt plink_log/all_gws_varids.txt > plink_log/flip.txt

# Flip alleles and check if the flipped alleles exist
while IFS=':' read -r -a varids; do
    new_line="${varids[0]}:${varids[1]}:${varids[3]}:${varids[2]}"
    echo "$new_line" >> plink_log/flipped_alleles.txt
done < plink_log/flip.txt

# Extract bim file for flipped variants of interest
plink2 \
--bfile /well/lindgren/samvida/Resources/1000Genomes/hg38/all_hg38 \
--allow-extra-chr \
--extract plink_log/flipped_alleles.txt \
--make-just-bim \
--threads 3 \
--memory 15000 \
--out plink_log/checking_flipped_alleles

# Collate list of variants that do exist in 1000G
awk '{print $2}' plink_log/checking_gws_varids.bim > plink_log/retain_1000G_variants.txt
awk '{print $2}' plink_log/checking_flipped_alleles.bim >> plink_log/retain_1000G_variants.txt

# Submit R script to get lead SNPs
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren/samvida/hormones_infertility/scripts/get_lead_snps_distance_pruning.R \
--hormone=${HORMONE} \
--sex_strata=${SEX_STRATA} \
--anc_group=${ANC_GROUP} 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
