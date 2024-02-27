#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J add_annotations_opentargets
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/add_annotations_opentargets-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/conditional_analysis

# Sort and get unique values from OpenTargets file
awk -v OFS="\t" 'NR > 1 {print "chr"$1":"$2, $3, $4, $5, $6}' \
/well/lindgren/samvida/Resources/OpenTargets/v2g_annots/v2g_annots_hg38.txt | sort -k1 | uniq > tmp_v2g_annots.txt

# Join with the classified lead SNPs file by position
join -t $'\t' -1 1 -2 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 2.2 2.3 2.4 2.5 \
<(sort -k1 classified_all_lead_snps_all_strata.txt) \
tmp_v2g_annots.txt \
> annotated_classified_all_lead_snps_all_strata.txt

# Add header
sed -i '1s/^/ID\tRSID\tCHROM\tGENPOS\tMAF\tAllele1\tAllele2\tFreq1\tFreqSE\tBETA\tSE\tPVALUE\tDirection\tHetPVal\thormone\tsex_strata\tancestry\tclassification\tVARID_1KG\tstrata\tgene_id\tfpred_max_label\tfpred_max_score\thgnc_symbol\n/' \
annotated_classified_all_lead_snps_all_strata.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

