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

export COMORMENT=/well/lindgren/samvida/Resources
export ANALYSIS_DIR=/well/lindgren/samvida/hormones_infertility
export SINGULARITY_BIND="$COMORMENT/mixer/reference:/REF:ro"
export SIF=$COMORMENT/mixer/singularity

singularity shell --bind ${ANALYSIS_DIR} --home $PWD:/home $SIF/mixer.sif

for pheno1 in female_infertility_analysis1 female_infertility_analysis3 female_infertility_analysis5; do
    for pheno2 in PCOS endometriosis hmb huber_twinning uterine_fibroids FSH_F Testosterone_F verdiesen_AMH williams_TSH Body_fat_percentage Body_mass_index Comparative_body_size_at_age_10 Waist_circumference Hip_circumference; do
        python /tools/mixer/precimed/mixer_figures.py combine --json results/$pheno1.fit.rep@.json  --out results/$pheno1.fit
        python /tools/mixer/precimed/mixer_figures.py combine --json results/$pheno1.test.rep@.json  --out results/$pheno1.test
        python /tools/mixer/precimed/mixer_figures.py combine --json results/$pheno2.fit.rep@.json  --out results/$pheno2.fit
        python /tools/mixer/precimed/mixer_figures.py combine --json results/$pheno2.test.rep@.json  --out results/$pheno2.test
        python /tools/mixer/precimed/mixer_figures.py combine --json results/${pheno1}_vs_${pheno2}.fit.rep@.json  --out results/${pheno1}_vs_${pheno2}.fit
        python /tools/mixer/precimed/mixer_figures.py combine --json results/${pheno1}_vs_${pheno2}.test.rep@.json  --out results/${pheno1}_vs_${pheno2}.test

        python /tools/mixer/precimed/mixer_figures.py one --json results/$pheno1.fit.json results/$pheno2.fit.json --out plots/${pheno1}_and_${pheno2}.fit --trait1 $pheno1 $pheno2 --statistic mean std --ext svg
        python /tools/mixer/precimed/mixer_figures.py one --json results/$pheno1.test.json results/$pheno2.test.json --out plots/${pheno1}_and_${pheno2}.test --trait1 $pheno1 $pheno2 --statistic mean std --ext svg

        python /tools/mixer/precimed/mixer_figures.py two --json-fit results/${pheno1}_vs_${pheno2}.fit.json --json-test results/${pheno1}_vs_${pheno2}.test.json --out plots/${pheno1}_vs_${pheno2} --trait1 $pheno1 --trait2 $pheno2 --statistic mean std --ext svg
    done  
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
