#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 16/03/23

# LD-scores from LDSC Bulik-Sullivan et al.

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J create_custom_baseline_ldscores
#SBATCH -o /well/lindgren/samvida/create_custom_baseline_ldscores-%j.out

echo "########################################################"
echo "Slurm Job ID: ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/Resources/LDSC

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

for CHR in {1..22}; do
	ldsc.py \
	--l2 \
	--bfile /well/lindgren/samvida/Resources/1000Genomes/chr_bedfiles_with_rsid_cm/1000G_Phase3_chr${CHR} \
	--ld-wind-cm 1 \
	--out baseline_ldscores/1000Genomes_Phase3_chr${CHR} \
	--print-snps /well/lindgren/samvida/Resources/HAPMAP3/hm3_chr${CHR}_only_snp.snplist
done

# cd /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility

# cat all_hormones_infertility_strata.txt | while read STRATA;
# do
# 	# European ancestry only, baseline
# 	ldsc.py \
# 	--h2 munged_sumstats/${STRATA}_EUR_for_ldsc.sumstats.gz \
# 	--ref-ld-chr /well/lindgren/samvida/Resources/LDSC/baseline_ldscores/1000Genomes_Phase3_chr \
# 	--w-ld-chr /well/lindgren/resources/sldsc/aux_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
# 	--overlap-annot \
# 	--frqfile-chr /well/lindgren/resources/sldsc/aux_files/1000G_Phase3_frq/1000G.EUR.QC. \
# 	--out partitioned_heritability/${STRATA}_EUR_baseline

# 	# All ancestries, baseline
# 	ldsc.py \
# 	--h2 munged_sumstats/${STRATA}_all_anc_for_ldsc.sumstats.gz \
# 	--ref-ld-chr /well/lindgren/samvida/Resources/LDSC/baseline_ldscores/1000Genomes_Phase3_chr \
# 	--w-ld-chr /well/lindgren/resources/sldsc/aux_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
# 	--overlap-annot \
# 	--frqfile-chr /well/lindgren/resources/sldsc/aux_files/1000G_Phase3_frq/1000G.EUR.QC. \
# 	--out partitioned_heritability/${STRATA}_all_anc_baseline
# done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
