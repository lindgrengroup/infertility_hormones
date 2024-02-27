#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 16/03/23

# LD-scores from LDSC Bulik-Sullivan et al.

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH --array 1-18
#SBATCH -J ovary_celltype_specific_partitioned_ldsc
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/ovary_celltype_specific_partitioned_ldsc-%j.out

echo "########################################################"
echo "Slurm Job ID: ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility
module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

STRATA=$( sed "${SLURM_ARRAY_TASK_ID}q;d" all_hormones_infertility_strata.txt )

# European ancestry only, partitioned LDSC heritability
ldsc.py \
--h2-cts munged_sumstats/${STRATA}_EUR_for_ldsc.sumstats.gz \
--ref-ld-chr /well/lindgren/samvida/Resources/LDSC/baseline_ldscores/1000Genomes_Phase3_chr \
--ref-ld-chr-cts /well/lindgren/samvida/hormones_infertility/ovary_celltype_markers_melody/for_ldsc/absolute_filepaths.ldcts \
--w-ld-chr /well/lindgren/resources/sldsc/aux_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--out ovary_celltypes_heritability/${STRATA}_EUR_ovary 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
