#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J get_reported_hormones_chrpos_hg38
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/get_reported_hormones_chrpos_hg38-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/Resources/GWASCatalog

# Extract variants associated with hormones
head -n 1 all_associations_230327.txt > gwascat_hormone_associations.txt
grep -F -w -f all_hormone_traits.txt all_associations_230327.txt >> gwascat_hormone_associations.txt

# Create bed file
Rscript scripts/make_gwascat_bed.R

# Get the variant IDs in 1KG hg38 format
module load PLINK/2.00a2.3_x86_64

plink2 \
--bfile /well/lindgren/samvida/Resources/1000Genomes/hg38/all_hg38 \
--allow-extra-chr \
--extract range gwascat_hormone_associations_hg38.bed \
--make-just-bim \
--threads 3 \
--memory 15000 \
--out gwascat_hormone_associations_hg38_1kg_variants

# Merge the 1KG variant ID into reported SNPs list

awk '{ print "chr"$1":"$4, $2 }' \
gwascat_hormone_associations_hg38_1kg_variants.bim | uniq | sort -k1 > gwascat_hormone_associations_hg38_1kg_variants.txt

cd /well/lindgren/samvida

# Get chromosome and position, sort by chr:pos
awk -F '\t' 'BEGIN {OFS = FS} { NR>1 ; print $1":"$3, $4}' \
Resources/GWASCatalog/gwascat_hormone_associations_hg38.bed | uniq | sort -k1 \
> hormones_infertility/conditional_analysis/published_snps/tmp.txt

join -1 1 -2 1 -o 1.1 2.2 1.2 \
hormones_infertility/conditional_analysis/published_snps/tmp.txt \
Resources/GWASCatalog/gwascat_hormone_associations_hg38_1kg_variants.txt \
> hormones_infertility/conditional_analysis/published_snps/all_hormones_reported.txt

# For each hormone, get a separate list of SNPs
declare -a STRATA=("FSH" "LH" "Oestradiol" "Progesterone" "Testosterone")
for strata_name in "${STRATA[@]}"; do
	grep -F -w -f hormones_infertility/conditional_analysis/published_snps/${strata_name}_traits.txt \
	hormones_infertility/conditional_analysis/published_snps/all_hormones_reported.txt \
	> hormones_infertility/conditional_analysis/published_snps/${strata_name}_reported.txt
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
