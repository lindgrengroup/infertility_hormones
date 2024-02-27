#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J extract_sig_hits
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/extract_sig_hits-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility
mkdir combine_wes_gwas/sig_hits

declare -a STRATA=("FSH_F" "FSH_M" "FSH_sex_comb" "LH_F" "LH_M" "LH_sex_comb" "Oestradiol_F" "Oestradiol_M" "Oestradiol_sex_comb" "Testosterone_F" "Testosterone_M" "Testosterone_sex_comb")

for strata_name in "${STRATA[@]}"; do
	# Extract GWS hits (P <= 5E-8) from common variant GWAS
	head -n 1 meta_results_230613/filtered/${strata_name}_EUR_filtered.txt > combine_wes_gwas/sig_hits/${strata_name}_gwas_sig_hits.txt
	awk -F'\t' 'BEGIN {OFS=FS} NR>1 { if ($11 <= 5E-8) print }' meta_results_230613/filtered/${strata_name}_EUR_filtered.txt \
	>> combine_wes_gwas/sig_hits/${strata_name}_gwas_sig_hits.txt
	# Get chr:pos alone for annots
	awk -F '\t' 'NR>1 { gsub(/chr23/, "chrX", $1); print $1}' combine_wes_gwas/sig_hits/${strata_name}_gwas_sig_hits.txt \
	> combine_wes_gwas/sig_hits/tmp_${strata_name}_gwas_chrpos.txt

	# Extract exome-wide significant variants (P <= 1E-7) from rare variant GWAS
	head -n 1 exome_seq_results/results_2310/cat_results/${strata_name}_variant.txt > combine_wes_gwas/sig_hits/${strata_name}_wes_variant_sig_hits.txt
	awk -F '\t' 'BEGIN {OFS=FS} NR>1 { if ($13 <= 1E-7) print }' exome_seq_results/results_2310/cat_results/${strata_name}_variant.txt \
	>> combine_wes_gwas/sig_hits/${strata_name}_wes_variant_sig_hits.txt
	# Get chr:pos alone for annots
	awk -F '\t' 'NR>1 { print "chr"$1":"$2":"}' combine_wes_gwas/sig_hits/${strata_name}_wes_variant_sig_hits.txt \
	> combine_wes_gwas/sig_hits/tmp_${strata_name}_wes_variant_chrpos.txt

	# Extract significant genes (P <= 5E-6) from (any) gene-based tests
	head -n 1 exome_seq_results/results_2310/cat_results/${strata_name}_gene.txt > combine_wes_gwas/sig_hits/${strata_name}_wes_gene_sig_hits.txt
	awk -F '\t' 'BEGIN {OFS=FS} NR>1 { if ($4 <= 5E-6) print }' exome_seq_results/results_2310/cat_results/${strata_name}_gene.txt \
	>> combine_wes_gwas/sig_hits/${strata_name}_wes_gene_sig_hits.txt
	awk -F '\t' 'NR>1 { print $1}' combine_wes_gwas/sig_hits/${strata_name}_wes_gene_sig_hits.txt \
	> combine_wes_gwas/sig_hits/tmp_${strata_name}_wes_gene.txt
done

# Extract annotation files for relevant variants
cd combine_wes_gwas/sig_hits

# Common variants
cat tmp_*_gwas_chrpos.txt | sort -u > all_hormones_gwas_chrpos.txt
head -n 1 /well/lindgren/samvida/Resources/OpenTargets/v2g_annots/v2g_annots_hg38.txt \
> all_hormones_sig_gwas_annots.txt 
grep -F -w -f all_hormones_gwas_chrpos.txt /well/lindgren/samvida/Resources/OpenTargets/v2g_annots/v2g_annots_hg19.txt \
>> all_hormones_sig_gwas_annots.txt 

# WES variants
cat tmp_*_wes_variant_chrpos.txt | sort -u > all_hormones_wes_variant_chrpos.txt
head -n 1 /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/ukb_wes_450k.july.qced.brava.v6.worst_csq_by_gene_canonical.tsv \
> all_hormones_sig_wes_variant_annots.txt 
grep -F -f all_hormones_wes_variant_chrpos.txt /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/ukb_wes_450k.july.qced.brava.v6.worst_csq_by_gene_canonical.tsv \
>> all_hormones_sig_wes_variant_annots.txt 

# Genes
cat tmp_*_wes_gene.txt | sort -u > all_hormones_wes_gene.txt
head -n 1 /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/gene_canonical.txt \
> all_hormones_sig_wes_gene_annots.txt 
grep -F -w -f all_hormones_wes_gene.txt /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/variant_annotations_for_wes/gene_canonical.txt \
>> all_hormones_sig_wes_gene_annots.txt 

rm tmp_*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
