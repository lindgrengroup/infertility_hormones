#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J lava_cross_trait_rg
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/all_phenos_lava_rg-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module purge

cd /well/lindgren/samvida/duncan_ldsc/ldsc
# conda env create --file environment.yml
source activate ldsc

cd /well/lindgren/samvida/hormones_infertility/lava_local

# Read filenames (all hormones, infertility, reproductive diseases, and obesity-associated traits)
mapfile -t sumstat_filenames < <(awk -F'\t' 'NR>1 {print $4}' input_sumstats.txt)

# Loop through all the files to make each one the target (one at a time), except the last file as it would have already been target
for ((target=0; target<${#sumstat_filenames[@]}-1; target++)); do
    # initialise
    all_phenos="${sumstat_filenames[$target]}"
    # loop through the remaining filenames
    for ((i=target+1; i<${#sumstat_filenames[@]}; i++)); do
        all_phenos="$all_phenos,${sumstat_filenames[$i]}"
    done

    /well/lindgren/samvida/duncan_ldsc/ldsc/ldsc.py \
    --rg $all_phenos \
    --ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
    --w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
    --out cross_trait_ldsc/${target}_vs_all_phenos_rg \
    --write-rg
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
