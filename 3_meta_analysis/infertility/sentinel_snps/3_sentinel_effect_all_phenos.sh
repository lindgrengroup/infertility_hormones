#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 07/08/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH --cpus-per-task=1
#SBATCH -J get_all_infertility_sentinels
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/get_all_infertility_sentinels-%j.out

echo `date`: Executing task ${SLURM_ARRAY_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/

# Compile a single list of sentinel SNPs for all male and female infertility
# female infertility
awk 'NR > 1 { print $1 }' infertility_meta_mvp/lead_snps/female_infertility_analysis1_all_sentinel_SNPs.txt \
> infertility_sumstats/tmp_female_infertility_all_sentinels.txt
for TRAITN in 2 3 4 5; do 
	awk 'NR > 1 { print $1 }' /well/lindgren/laura/projects/infertility/sentinel_snps/female_infertility_analysis${TRAITN}_all_sentinel_SNPs.txt \
	>> infertility_sumstats/tmp_female_infertility_all_sentinels.txt
done
sort infertility_sumstats/tmp_female_infertility_all_sentinels.txt | uniq > infertility_sumstats/female_infertility_all_sentinels.txt
rm infertility_sumstats/tmp_*

# male infertility
awk 'NR > 1 { print $1 }' infertility_meta_mvp/lead_snps/male_infertility_all_sentinel_SNPs.txt \
> infertility_sumstats/male_infertility_all_sentinels.txt

# Extract summary statistics for all sentinels
head -n 1 infertility_meta_mvp/filtered/female_infertility_analysis1_all_filtered.txt > infertility_sumstats/all_lead_snps_female_infertility_analysis1_all.txt
grep -Fwf infertility_sumstats/female_infertility_all_sentinels.txt infertility_meta_mvp/filtered/female_infertility_analysis1_all_filtered.txt \
>> infertility_sumstats/all_lead_snps_female_infertility_analysis1_all.txt

head -n 1 /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/female_infertility_analysis1_eur_MA_results_chr_pos.txt > infertility_sumstats/all_lead_snps_female_infertility_analysis1_eur.txt
grep -Fwf infertility_sumstats/female_infertility_all_sentinels.txt /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/female_infertility_analysis1_eur_MA_results_chr_pos.txt \
>> infertility_sumstats/all_lead_snps_female_infertility_analysis1_eur.txt

for TRAITN in 2 3 4 5; do
	for anc in all eur; do
		head -n 1 /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/female_infertility_analysis${TRAITN}_${anc}_MA_results_chr_pos.txt \
		> infertility_sumstats/all_lead_snps_female_infertility_analysis${TRAITN}_${anc}.txt
		grep -Fwf infertility_sumstats/female_infertility_all_sentinels.txt \
		/well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/female_infertility_analysis${TRAITN}_${anc}_MA_results_chr_pos.txt \
		>> infertility_sumstats/all_lead_snps_female_infertility_analysis${TRAITN}_${anc}.txt
	done
done

head -n 1 infertility_meta_mvp/filtered/male_infertility_all_filtered.txt > infertility_sumstats/all_lead_snps_male_infertility_all.txt
grep -Fwf infertility_sumstats/male_infertility_all_sentinels.txt infertility_meta_mvp/filtered/male_infertility_all_filtered.txt \
>> infertility_sumstats/all_lead_snps_male_infertility_all.txt

head -n 1 /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/male_infertility_eur_MA_results_chr_pos.txt \
> infertility_sumstats/all_lead_snps_male_infertility_eur.txt
grep -Fwf infertility_sumstats/male_infertility_all_sentinels.txt /well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process/male_infertility_eur_MA_results_chr_pos.txt \
>> infertility_sumstats/all_lead_snps_male_infertility_eur.txt

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript scripts/plot_forest_sentinel_snps.R

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

