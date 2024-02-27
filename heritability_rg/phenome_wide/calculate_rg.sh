#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J calculate_rg
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/all_phenos_rg-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module purge

cd /well/lindgren/samvida/duncan_ldsc/ldsc
# conda env create --file environment.yml
source activate ldsc

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

cd /well/lindgren/samvida/hormones_infertility/rg_phenome_wide

sumstat_filenames=(/well/lindgren/resources/Neale_LDSC_all_phenos/*.tsv.bgz)
with_target_filenames=("/well/lindgren/samvida/hormones_infertility/rg_hormones_infertility/munged_sumstats/${STRATA}_EUR_for_ldsc.sumstats.gz" "${sumstat_filenames[@]}")
all_phenos=$(IFS=, ; echo "${with_target_filenames[*]}")

/well/lindgren/samvida/duncan_ldsc/ldsc/ldsc.py \
--rg $all_phenos \
--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
--out ${STRATA}_vs_all_phenos \
--write-rg

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
