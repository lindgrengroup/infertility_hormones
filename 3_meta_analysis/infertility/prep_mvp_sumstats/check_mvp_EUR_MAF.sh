#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 02/08/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J check_mvp_eur_maf
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/check_mvp_eur_maf-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_ARRAY_TASK_ID" 
echo "Slurm Task ID: $SLURM_ARRAY_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

# Join the old meta-analysis and MVP sumstats by ID

cd /well/lindgren/samvida/hormones_infertility

# Join by chr:pos
join -t $'\t' -1 1 -2 1 -o 1.1,1.7,2.11 \
<(sort -k1 public_sumstats/female_infertility_analysis1/filtered/MVP_female_infertility_analysis1_EUR_filtered.txt) \
<(sort -k1 /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/female_infertility_analysis1_eur_MA_results_chr_pos.txt) \
> infertility_meta_mvp/debug/compare_mafs_laura_mvp.txt

# Join this with the MA with MVP to see if the FreqSE is higher because of inclusion of MVP
join -t $'\t' -1 1 -2 1 -o 1.1,1.2,1.3,2.8 \
infertility_meta_mvp/debug/compare_mafs_laura_mvp.txt \
<(sort -k1 infertility_meta_mvp/filtered/female_infertility_analysis1_EUR_filtered.txt) \
> infertility_meta_mvp/debug/compare_mafs_laura_mvp_full.txt

echo -e "ID\tMAF_MVP\tMAF_old_meta\tFreqSE_new_meta" | cat - infertility_meta_mvp/debug/compare_mafs_laura_mvp_full.txt > temp && mv temp infertility_meta_mvp/debug/compare_mafs_laura_mvp_full.txt

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript scripts/check_mvp_EUR_MAF.R

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
