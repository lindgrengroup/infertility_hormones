#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J compare_public_no_public
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/compare_public_no_public-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/meta_results_230613_no_public

cat ../with_public_strata_list.txt | while read STRATA; do
	# Get lead SNPs list
	awk 'NR > 1 {print $1}' ../meta_results_230613/lead_snps/all_lead_snp_sumstats_${STRATA}_with_rsids.txt > \
	compare_public/${STRATA}_lead_snps.txt
	
	sed -i 's/chr23:/chrX:/g' ${STRATA}_1.tbl 
	grep -Fwf compare_public/${STRATA}_lead_snps.txt ${STRATA}_1.tbl > \
	compare_public/lead_snps_no_public_${STRATA}.txt

	sed 's/chr23:/chrX:/g' ../meta_results_230613/${STRATA}_1.tbl > ../meta_results_230613/tmp_${STRATA}.tbl 
	grep -Fwf compare_public/${STRATA}_lead_snps.txt ../meta_results_230613/tmp_${STRATA}.tbl > \
	compare_public/lead_snps_with_public_${STRATA}.txt

	rm ../meta_results_230613/tmp_${STRATA}.tbl 
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
