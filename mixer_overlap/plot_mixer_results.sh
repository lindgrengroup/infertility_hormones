#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 01/08/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J plot_mixer_results
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/plot_mixer_results-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/mixer_overlap

export ANALYSIS_DIR=/well/lindgren/samvida/hormones_infertility
singularity shell --bind ${ANALYSIS_DIR} --home $PWD:/home $SIF/mixer.sif

python /tools/mixer/precimed/mixer_figures.py combine --json SCZ.fit.rep@.json  --out SCZ.fit
python /tools/mixer/precimed/mixer_figures.py combine --json SCZ.test.rep@.json  --out SCZ.test
python /tools/mixer/precimed/mixer_figures.py combine --json INT.fit.rep@.json  --out INT.fit
python /tools/mixer/precimed/mixer_figures.py combine --json INT.test.rep@.json  --out INT.test
python /tools/mixer/precimed/mixer_figures.py combine --json SCZ_vs_INT.fit.rep@.json  --out SCZ_vs_INT.fit
python /tools/mixer/precimed/mixer_figures.py combine --json SCZ_vs_INT.test.rep@.json  --out SCZ_vs_INT.test

python /tools/mixer/precimed/mixer_figures.py one --json SCZ.fit.json INT.fit.json --out SCZ_and_INT.fit --trait1 SCZ INT --statistic mean std --ext svg
python /tools/mixer/precimed/mixer_figures.py one --json SCZ.test.json INT.test.json --out SCZ_and_INT.test --trait1 SCZ INT --statistic mean std --ext svg

python /tools/mixer/precimed/mixer_figures.py two --json-fit SCZ_vs_INT.fit.json --json-test SCZ_vs_INT.test.json --out SCZ_vs_INT --trait1 SCZ --trait2 INT --statistic mean std --ext svg

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
