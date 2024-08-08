#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 05/08/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J extract_sig_hits
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/extract_sig_hits-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility

# there are only three strata with significant WES variants
declare -a STRATA=("Oestradiol_F" "Testosterone_F" "Testosterone_M")

for strata_name in "${STRATA[@]}"; do
	# Get the sig variants in WES (just chr:pos)
	awk -F '\t' 'NR>1 { print "chr"$1":"$2}' combine_wes_gwas/sig_hits/${strata_name}_wes_variant_sig_hits.txt \
	> combine_wes_gwas/sig_hits/tmp_${strata_name}_wes_variant_chrpos.txt

	# Extract these hits from common variant GWAS
	head -n 1 meta_results_230613_no_ukb/filtered/${strata_name}_EUR_filtered.txt \
	> combine_wes_gwas/sig_hits/${strata_name}_gwas_sumstats_wes_variants.txt
	grep -wFf combine_wes_gwas/sig_hits/tmp_${strata_name}_wes_variant_chrpos.txt \
	meta_results_230613_no_ukb/filtered/${strata_name}_EUR_filtered.txt >> combine_wes_gwas/sig_hits/${strata_name}_gwas_sumstats_wes_variants.txt
done

rm combine_wes_gwas/sig_hits/tmp_*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
