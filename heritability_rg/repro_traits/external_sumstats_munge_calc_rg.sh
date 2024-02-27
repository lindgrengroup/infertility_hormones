#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J munge_sumstats_ldsc
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/repro_external_traits_munge_sumstats_ldsc-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/rg_repro_traits

# For twinning
# Create a chr:pos column, make the alleles capital letters, 
# subset to only the columns that are needed for LDSC, and sort by the chr:pos column
awk -v OFS="\t" '{print "chr"$3, toupper($4), toupper($5), $6, $7, $9}' \
external_sumstats/huber_twinning_010620.dat \
| sort -k1 > munged_sumstats/tmp_huber_twinning_sumstats.txt

join -t $'\t' -1 1 -2 6 -o 2.3,1.2,1.3,1.4,1.5,1.6,1.7,2.4,2.5 \
munged_sumstats/tmp_huber_twinning_sumstats.txt \
<(sort -k6 /well/lindgren/samvida/Resources/hg19/grch37.withchrpos.hm3.rsid.txt) \
> munged_sumstats/tmp_huber_twinning_sumstats_with_rsids.txt

# Remove rows where the map alleles don't match the sumstat alleles
awk 'BEGIN { FS = OFS = "\t" } { split($9, values, ","); for (i in values) if ($2 == $8 && $3 == values[i] || $3 == $8 && $2 == values[i]) { print $1,$2,$3,$4,$5,$6,$7; break } }' \
munged_sumstats/tmp_huber_twinning_sumstats_with_rsids.txt \
> munged_sumstats/huber_twinning_sumstats_with_rsids.txt

# Add header column
sed -i '1s/^/RSID\tAllele1\tAllele2\tFreq1\tBETA\tPVALUE\n/' \
munged_sumstats/huber_twinning_sumstats_with_rsids.txt
rm munged_sumstats/tmp_huber_twinning_*.txt

# For AMH
# Create a chr:pos column, make the alleles letters removing quotes, 
# subset to only the columns that are needed for LDSC, and sort by the chr:pos column
awk -v OFS="\t" '{gsub(/"/, "", $4); gsub(/"/, "", $5); print "chr"$1":"$2, $4, $5, $6, $7, $3}' \
external_sumstats/verdiesen_AMH_231127_buildGRCh37.tsv \
| sort -k1 > munged_sumstats/tmp_verdiesen_AMH_sumstats.txt

join -t $'\t' -1 1 -2 6 -o 2.3,1.2,1.3,1.4,1.5,1.6,1.7,2.4,2.5 \
munged_sumstats/tmp_verdiesen_AMH_sumstats.txt \
<(sort -k6 /well/lindgren/samvida/Resources/hg19/grch37.withchrpos.hm3.rsid.txt) \
> munged_sumstats/tmp_verdiesen_AMH_sumstats_with_rsids.txt

# Remove rows where the map alleles don't match the sumstat alleles
awk 'BEGIN { FS = OFS = "\t" } { split($9, values, ","); for (i in values) if ($2 == $8 && $3 == values[i] || $3 == $8 && $2 == values[i]) { print $1,$2,$3,$4,$5,$6,$7; break } }' \
munged_sumstats/tmp_verdiesen_AMH_sumstats_with_rsids.txt \
> munged_sumstats/verdiesen_AMH_sumstats_with_rsids.txt

# Add header column
sed -i '1s/^/RSID\tAllele1\tAllele2\tFreq1\tBETA\tPVALUE\n/' \
munged_sumstats/verdiesen_AMH_sumstats_with_rsids.txt
rm munged_sumstats/tmp_verdiesen_AMH_*.txt

# For TSH
# Create a chr:pos column, 
# subset to only the columns that are needed for LDSC, and sort by the chr:pos column
awk -v OFS="\t" '{print "chr"$1":"$2, $3, $4, $7, $5, $8}' \
external_sumstats/williams_TSH_231127.tsv \
| sort -k1 > munged_sumstats/tmp_williams_TSH_sumstats.txt

join -t $'\t' -1 1 -2 6 -o 2.3,1.2,1.3,1.4,1.5,1.6,1.7,2.4,2.5 \
munged_sumstats/tmp_williams_TSH_sumstats.txt \
<(sort -k6 /well/lindgren/samvida/Resources/hg19/grch37.withchrpos.hm3.rsid.txt) \
> munged_sumstats/tmp_williams_TSH_sumstats_with_rsids.txt

# Remove rows where the map alleles don't match the sumstat alleles
awk 'BEGIN { FS = OFS = "\t" } { split($9, values, ","); for (i in values) if ($2 == $8 && $3 == values[i] || $3 == $8 && $2 == values[i]) { print $1,$2,$3,$4,$5,$6,$7; break } }' \
munged_sumstats/tmp_williams_TSH_sumstats_with_rsids.txt \
> munged_sumstats/williams_TSH_sumstats_with_rsids.txt

# Add header column
sed -i '1s/^/RSID\tAllele1\tAllele2\tFreq1\tBETA\tPVALUE\n/' \
munged_sumstats/williams_TSH_sumstats_with_rsids.txt
rm munged_sumstats/tmp_williams_TSH_*.txt

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

# For twinning
munge_sumstats.py \
--sumstats munged_sumstats/huber_twinning_sumstats_with_rsids.txt \
--N-cas 8265 --N-con 264567 \
--snp RSID --a1 ALLELE1 --a2 ALLELE2 --p PVALUE \
--frq FREQ1 --signed-sumstats BETA,0 \
--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
--out /well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/huber_twinning_for_ldsc

# For AMH
munge_sumstats.py \
--sumstats munged_sumstats/verdiesen_AMH_sumstats_with_rsids.txt \
--N 7049 \
--snp RSID --a1 ALLELE1 --a2 ALLELE2 --p PVALUE \
--frq FREQ1 --signed-sumstats BETA,0 \
--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
--out /well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/verdiesen_AMH_for_ldsc

# For TSH
munge_sumstats.py \
--sumstats munged_sumstats/williams_TSH_sumstats_with_rsids.txt \
--N 247107 \
--snp RSID --a1 ALLELE1 --a2 ALLELE2 --p PVALUE \
--frq FREQ1 --signed-sumstats BETA,0 \
--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
--out /well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/williams_TSH_for_ldsc

module purge

cd /well/lindgren/samvida/duncan_ldsc/ldsc
# conda env create --file environment.yml
source activate ldsc

# Run through each infertility to calculate rG
for INFERT_STRATA in female_infertility_analysis1 female_infertility_analysis2 female_infertility_analysis3 female_infertility_analysis4 female_infertility_analysis5; do
	/well/lindgren/samvida/duncan_ldsc/ldsc/ldsc.py \
	--rg /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility/munged_sumstats/${INFERT_STRATA}_EUR_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/endometriosis_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/hmb_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/PCOS_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/uterine_fibroids_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/huber_twinning_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/verdiesen_AMH_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/williams_TSH_for_ldsc.sumstats.gz \
	--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--out /well/lindgren/samvida/hormones_infertility/rg_repro_traits/${INFERT_STRATA}_vs_repro_phenos \
	--write-rg
done

# Run through each hormone to calculate rG
for HORMONE_STRATA in FSH LH Oestradiol Testosterone; do
	/well/lindgren/samvida/duncan_ldsc/ldsc/ldsc.py \
	--rg /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility/munged_sumstats/${HORMONE_STRATA}_F_EUR_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/endometriosis_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/hmb_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/PCOS_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/uterine_fibroids_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/huber_twinning_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/verdiesen_AMH_for_ldsc.sumstats.gz,/well/lindgren/samvida/hormones_infertility/rg_repro_traits/munged_sumstats/williams_TSH_for_ldsc.sumstats.gz \
	--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--out /well/lindgren/samvida/hormones_infertility/rg_repro_traits/${HORMONE_STRATA}_vs_repro_phenos \
	--write-rg
done


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0



