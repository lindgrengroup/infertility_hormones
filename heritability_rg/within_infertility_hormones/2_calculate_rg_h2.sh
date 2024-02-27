#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J calculate_rg
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/rg_between_traits-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility
module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

# Strata to loop through
readarray -t STRATA < all_hormones_infertility_strata.txt
nstrat=${#STRATA[@]}
nstrat_less=$(expr $nstrat - 1)

for i in $(seq 1 $nstrat); do
	for j in $(seq $i $nstrat_less); do
		# European ancestry only
		ldsc.py \
		--rg munged_sumstats/${STRATA[$i-1]}_EUR_for_ldsc.sumstats.gz,munged_sumstats/${STRATA[$j]}_EUR_for_ldsc.sumstats.gz \
		--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--out results/rg_EUR_${STRATA[$i-1]}_${STRATA[$j]}

		# All ancestries
		ldsc.py \
		--rg munged_sumstats/${STRATA[$i-1]}_all_anc_for_ldsc.sumstats.gz,munged_sumstats/${STRATA[$j]}_all_anc_for_ldsc.sumstats.gz \
		--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--out results/rg_all_anc_${STRATA[$i-1]}_${STRATA[$j]}
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

