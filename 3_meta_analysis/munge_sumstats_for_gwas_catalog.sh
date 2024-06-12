#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 02/05/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 2
#SBATCH -J munge_sumstats_gwascat
#SBATCH -o /well/lindgren/samvida/hormones_infertility/sumstats_gwas_catalog/logs/munge-sumstats-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load Python/3.10.8-GCCcore-12.2.0
source /well/lindgren/samvida/python/gwas-catalog-sumstats-${MODULE_CPU_TYPE}/bin/activate

cd /well/lindgren/samvida/hormones_infertility

# INFERTILITY

for anc in all eur; do
# FEMALE INFERTILITY
	for analysis_strata in {1..5}; do
		# old file: 1.ID(chrN:pos) 2.RSID 3.Allele1 4.Allele2 5.Freq1 6.FreqSE 7.MinFreq 8.MaxFreq 9.BETA 10.SE 11.PVALUE 12.Direction 13.HetISq 14.HetChiSq 15.HetDf 16.HetPVal 17.N_CASES 18.N_CONTROLS
		# new file: chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value variant_id rs_id n_cases n_controls direction HetISq HetChiSq HetDf HetPVal
		awk 'BEGIN{FS=" "; OFS="\t"} {
			if ($1 ~ /^chr[0-9]+:[0-9]+$/) {
				# Extract chrN and [position] from the first column
				match($1, /^chr([0-9]+):([0-9]+)$/, arr)
				CHR = arr[1]
				POS = arr[2]

				# Print the extracted columns and the remaining columns in a new arrangement
				print CHR, POS, $3, $4, $9, $10, $5, $11, $2, $17, $18, $12, $13, $14, $15, $16 
			}
		}' infertility_sumstats/female_infertility_analysis${analysis_strata}_${anc}_with_rsids.txt \
		> sumstats_gwas_catalog/female_infertility_analysis${analysis_strata}_${anc}_sumstats.tsv

		sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id\tn_cases\tn_controls\tdirection\tHetISq\tHetChiSq\tHetDf\tHetPVal/' \
		sumstats_gwas_catalog/female_infertility_analysis${analysis_strata}_${anc}_sumstats.tsv

		# Validate
		gwas-ssf validate --errors-out sumstats_gwas_catalog/female_infertility_analysis${analysis_strata}_${anc}_sumstats.tsv
		# gzip
		gzip sumstats_gwas_catalog/female_infertility_analysis${analysis_strata}_${anc}_sumstats.tsv
	done

# MALE INFERTILITY
	awk 'BEGIN{FS=" "; OFS="\t"} {
		if ($1 ~ /^chr[0-9]+:[0-9]+$/) {
			# Extract chrN and [position] from the first column
			match($1, /^chr([0-9]+):([0-9]+)$/, arr)
			CHR = arr[1]
			POS = arr[2]

			# Print the extracted columns and the remaining columns in a new arrangement
			print CHR, POS, $3, $4, $9, $10, $5, $11, $2, $17, $18, $12, $13, $14, $15, $16 
		}
	}' infertility_sumstats/male_infertility_${anc}_with_rsids.txt \
	> sumstats_gwas_catalog/male_infertility_${anc}_sumstats.tsv
	
	sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id\tn_cases\tn_controls\tdirection\tHetISq\tHetChiSq\tHetDf\tHetPVal/' \
	sumstats_gwas_catalog/male_infertility_${anc}_sumstats.tsv

	# Validate
	gwas-ssf validate --errors-out sumstats_gwas_catalog/male_infertility_${anc}_sumstats.tsv
	# gzip
	gzip sumstats_gwas_catalog/male_infertility_${anc}_sumstats.tsv
done

# HORMONES

for HR in FSH LH Oestradiol Progesterone Testosterone; do
	for SS in F M sex_comb; do
		for ANC in EUR all_anc; do
			# old file: 1.ID 2.RSID 3.CHROM 4.GENPOS 5.MAF 6.Allele1 7.Allele2 8.Freq1 9.FreqSE 10.BETA 11.SE 12.PVALUE 13.Direction 14.HetPVal
			# new file: chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rs_id direction HetPVal
			awk 'BEGIN{FS=OFS="\t"} { 
				print $3, $4, $6, $7, $10, $11, $8, $12, $2, $13, $14
			}' /well/lindgren/samvida/hormones_infertility/meta_results_230613/add_rsids/${HR}_${SS}_${ANC}_with_rsids.txt \
			> sumstats_gwas_catalog/${HR}_${SS}_${ANC}_sumstats.tsv
			
			sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id\tdirection\tHetPVal/' \
			sumstats_gwas_catalog/${HR}_${SS}_${ANC}_sumstats.tsv

			# Validate
			gwas-ssf validate --errors-out sumstats_gwas_catalog/${HR}_${SS}_${ANC}_sumstats.tsv
			# gzip
			gzip sumstats_gwas_catalog/${HR}_${SS}_${ANC}_sumstats.tsv
		done
	done
done

# OTHER REPRODUCTIVE DISEASES META-ANALYSED (PCOS, heavy menstrual bleeding, endometriosis, uterine fibroids)

for REPRO in PCOS hmb endometriosis uterine_fibroids; do
	# old file: 1.RSID 2.CHR 3.POS 4.ALLELE1 5.ALLELE2 6.FREQ1 7.BETA 8.SE 9.PVALUE 10.NCASE 11.NCONTROL 
	# new file: chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rs_id
	awk 'BEGIN{FS=OFS="\t"} { 
		# Replace chrX with chr23
		gsub(/^chrX/, "chr23", $2)
		# Strip "chr"
		sub(/^chr/, "", $2)
		# Only print rows with numeric chromosomes
		if ($2 ~ /^[0-9]+$/) {
			print $2, $3, $4, $5, $7, $8, $6, $9, $1
		}
	}' /well/lindgren/samvida/hormones_infertility/colocalisation/sentinel_windows/liftover/${REPRO}_hg38_with_samplesizes.txt \
	> sumstats_gwas_catalog/${REPRO}_sumstats.tsv
	
	sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id/' \
	sumstats_gwas_catalog/${REPRO}_sumstats.tsv

	# Validate
	gwas-ssf validate --errors-out sumstats_gwas_catalog/${REPRO}_sumstats.tsv
	# gzip
	gzip sumstats_gwas_catalog/${REPRO}_sumstats.tsv
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
