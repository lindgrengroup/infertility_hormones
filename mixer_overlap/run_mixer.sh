#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 01/08/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH --cpus-per-task=8
#SBATCH --array=1-20
#SBATCH -J run_mixer
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/run_mixer-%j.out

echo `date`: Executing task ${SLURM_ARRAY_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/mixer_overlap

# pass variables for: pheno1, pheno2, pheno1-sumstats-location, pheno2-sumstats-location
echo "passing variables..."
covars=${covars//|/,}
echo ${covars}

export COMORMENT=/well/lindgren/samvida/Resources
export ANALYSIS_DIR=/well/lindgren/samvida/hormones_infertility
export SINGULARITY_BIND="$COMORMENT/mixer/reference:/REF:ro"
export SIF=$COMORMENT/mixer/singularity
export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 8"
export REP="rep${SLURM_ARRAY_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps"

export PYTHON="singularity exec --bind ${ANALYSIS_DIR} --home=$PWD:/home $SIF/mixer.sif python"

$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file ${pheno1_loc} --out results/${pheno1}.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file ${pheno2_loc} --out results/${pheno2}.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file ${pheno1_loc} --trait2-file ${pheno2_loc} --trait1-params results/${pheno1}.fit.$REP.json --trait2-params results/${pheno2}.fit.$REP.json --out results/${pheno1}_vs_${pheno2}.fit.$REP

$PYTHON /tools/mixer/precimed/mixer.py test1 $MIXER_COMMON_ARGS --trait1-file ${pheno1_loc} --load-params results/${pheno1}.fit.$REP.json --out results/${pheno1}.test.$REP
$PYTHON /tools/mixer/precimed/mixer.py test1 $MIXER_COMMON_ARGS --trait1-file ${pheno2_loc} --load-params results/${pheno2}.fit.$REP.json --out results/${pheno2}.test.$REP
$PYTHON /tools/mixer/precimed/mixer.py test2 $MIXER_COMMON_ARGS --trait1-file ${pheno1_loc} --trait2-file ${pheno2_loc} --load-params results/${pheno1}_vs_${pheno2}.fit.$REP.json --out results/${pheno1}_vs_${pheno2}.test.$REP

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
