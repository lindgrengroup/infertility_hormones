#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J munge_sumstats_ldsc
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/repro_traits_munge_sumstats_ldsc-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

# For repro_traits
while IFS=$'\t' read -r -a STRATA_INFO
do
	munge_sumstats.py \
	--sumstats /well/lindgren-ukbb/projects/ukbb-11867/samvida/obesity_wrh_thesis/two_sample_mr/repro_outcomes/${STRATA_INFO[0]}_sumstats_with_rsids.txt \
	--N-cas ${STRATA_INFO[1]} --N-con ${STRATA_INFO[2]} \
	--snp RSID --a1 ALLELE1 --a2 ALLELE2 --p PVALUE \
	--frq FREQ1 --signed-sumstats BETA,0 \
	--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
	--out /well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/${STRATA_INFO[0]}_for_ldsc
done < /well/lindgren/samvida/hormones_infertility/rg_repro_traits/repro_files_list.txt

module purge

cd /well/lindgren/samvida/duncan_ldsc/ldsc
# conda env create --file environment.yml
source activate ldsc

# Run through each infertility to calculate rG
for INFERT_STRATA in female_infertility_analysis1 female_infertility_analysis2 female_infertility_analysis3 female_infertility_analysis4 female_infertility_analysis5; do
	/well/lindgren/samvida/duncan_ldsc/ldsc/ldsc.py \
	--rg /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility/munged_sumstats/${INFERT_STRATA}_EUR_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/endometriosis_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/hmb_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/PCOS_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/uterine_fibroids_for_ldsc.sumstats.gz \
	--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--out /well/lindgren/samvida/hormones_infertility/rg_repro_traits/${INFERT_STRATA}_vs_repro_phenos \
	--write-rg
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0



