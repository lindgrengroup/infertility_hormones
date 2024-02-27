#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 31/10/22

# LiftOver all sumstats from hg19 to hg38

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J liftover_all_hg19_to_hg38
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/liftover_all_hg19_to_hg38-%j.out

echo "########################################################"
echo "Slurm Job ID: ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load liftOver/20210519

# ALSPAC
cd /well/lindgren/samvida/hormones_infertility/data/ALSPAC
mkdir lifted
mkdir lifted/tmp_files
touch lifted/liftover_hg19Tohg38.log
cat pheno_list.txt | while read PHENO;
do
	if [ -f "alsp_${PHENO}_gwas_res_ae.txt" ]; then
		printf "** Phenotype: ${PHENO} in sex strata: females\n" >> lifted/liftover_hg19Tohg38.log

		# Replace 23 with X for chrX
		awk -v FS="\t" -v OFS="\t" '{ sub("23", "X", $2) } 1' alsp_${PHENO}_gwas_res_ae.txt > tmp.txt
		# Add unique_identifier to sumstats (chr_pos_A1_A2)
		awk -v FS="\t" -v OFS="\t" '{ print "chr"$2"_"$3"_"$5"_"$6, $0 }' tmp.txt \
		> with_unique_id_alsp_${PHENO}_gwas_res_ae.txt
		rm tmp.txt

		# Create bed file with chrX, pos0, pos1, unique_id (CHR_POS_A1_A2)
		awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$3, $4-1, $4, $1}' \
		with_unique_id_alsp_${PHENO}_gwas_res_ae.txt \
		> ${PHENO}_hg19.bed
			
		liftOver ${PHENO}_hg19.bed \
		/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
		${PHENO}_hg38.bed \
		${PHENO}_unlifted.bed

		# Log the % of SNPs succesfully lifted over
		hg19=$(cat ${PHENO}_hg19.bed | wc -l)
		printf "\t # SNPs in hg19: $hg19 \n" >> lifted/liftover_hg19Tohg38.log
		hg38=$(cat ${PHENO}_hg38.bed | wc -l)
		printf "\t # SNPs lifted over to hg38: $hg38 \n" >> lifted/liftover_hg19Tohg38.log
		perc_lifted=$(echo "scale=3;${hg38}/${hg19}" | bc)
		printf "\t percent successfully lifted: $(echo "scale=3;${hg38}/${hg19}*100" | bc) \n" >> lifted/liftover_hg19Tohg38.log

		# Create a new sumstats file by joining via the unique_id
		# unique_id: 4th column in the lifted over bed file and 1st column in the original sumstats
		# Add header
		head -1 with_unique_id_alsp_${PHENO}_gwas_res_ae.txt \
		> lifted/${PHENO}_hg38.txt
		join -1 4 -2 1 -o 2.1,2.2,1.1,1.3,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15 \
		<(sort -k 4 ${PHENO}_hg38.bed) <(sort -k 1 with_unique_id_alsp_${PHENO}_gwas_res_ae.txt) \
		> lifted/tmp_files/alsp_${PHENO}_tmp_joined.txt
		awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15 }' \
		lifted/tmp_files/alsp_${PHENO}_tmp_joined.txt \
		>> lifted/${PHENO}_hg38.txt
		mv *.bed lifted/tmp_files
	fi
done 

# EstBB
cd /well/lindgren/samvida/hormones_infertility/data/EstBB
mkdir lifted
mkdir lifted/tmp_files
touch lifted/liftover_hg19Tohg38.log
cat pheno_list.txt | while read PHENO;
do
	for ss in female male combined; do
		if [ -f "${PHENO}_${ss}_EstBB_European_JF_05012023.txt" ]; then
			printf "** Phenotype: ${PHENO} in sex strata: ${ss}\n" >> lifted/liftover_hg19Tohg38.log

			# Replace 23 with X for chrX
			awk -v FS=" " -v OFS=" " '{ sub("23", "X", $1) } 1' ${PHENO}_${ss}_EstBB_European_JF_05012023.txt > tmp.txt
			# Add unique_identifier to sumstats (chr_pos_A1_A2)
			awk -v FS=" " -v OFS=" " '{ print "chr"$1"_"$2"_"$4"_"$5, $0 }' tmp.txt \
			> with_unique_id_${PHENO}_${ss}_EstBB_European_JF_05012023.txt
			rm tmp.txt

			# Create bed file with chrX, pos0, pos1, unique_id
			awk -v FS=" " -v OFS="\t" 'NR > 1 {print "chr"$2, $3-1, $3, $1}' \
			with_unique_id_${PHENO}_${ss}_EstBB_European_JF_05012023.txt \
			> ${PHENO}_${ss}_hg19.bed
				
			liftOver ${PHENO}_${ss}_hg19.bed \
			/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
			${PHENO}_${ss}_hg38.bed \
			${PHENO}_${ss}_unlifted.bed

			# Log the % of SNPs succesfully lifted over
			hg19=$(cat ${PHENO}_${ss}_hg19.bed | wc -l)
			printf "\t # SNPs in hg19: $hg19 \n" >> lifted/liftover_hg19Tohg38.log
			hg38=$(cat ${PHENO}_${ss}_hg38.bed | wc -l)
			printf "\t # SNPs lifted over to hg38: $hg38 \n" >> lifted/liftover_hg19Tohg38.log
			perc_lifted=$(echo "scale=3;${hg38}/${hg19}" | bc)
			printf "\t percent successfully lifted: $(echo "scale=3;${hg38}/${hg19}*100" | bc) \n" >> lifted/liftover_hg19Tohg38.log

			# Create a new sumstats file by joining via unique_id
			# unique_id: 4th column in the lifted over bed file and 1st column in the original sumstats
			# Add header
			head -1 with_unique_id_${PHENO}_${ss}_EstBB_European_JF_05012023.txt \
			> lifted/${PHENO}_${ss}_hg38.txt
			join -1 4 -2 1 -o 2.1,1.1,1.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15 \
			<(sort -k 4 ${PHENO}_${ss}_hg38.bed) <(sort -k 1 with_unique_id_${PHENO}_${ss}_EstBB_European_JF_05012023.txt) \
			>> lifted/${PHENO}_${ss}_hg38.txt
			mv *.bed lifted/tmp_files
		fi
	done
done 

# UKBB
cd /well/lindgren/samvida/hormones_infertility/data/UKBB
mkdir lifted
mkdir lifted/tmp_files
touch lifted/liftover_hg19Tohg38.log
cat pheno_list.txt | while read PHENO;
do
	for anc in non_finnish_EUR AFR EAS SAS; do
		for ss in F M sex_comb; do
			if [ -f "filtered_hg19/${PHENO}_${anc}_${ss}_filtered.txt" ]; then
				printf "** Phenotype: ${PHENO} in ancestry ${anc} and sex strata ${ss}\n" >> lifted/liftover_hg19Tohg38.log

				# Replace 23 with X for chrX
				awk -v FS="\t" -v OFS="\t" '{ sub("23", "X", $2) } 1' filtered_hg19/${PHENO}_${anc}_${ss}_filtered.txt > tmp.txt
				# Add unique_identifier to sumstats (chr_pos_A1_A2)
				awk -v FS="\t" -v OFS="\t" '{ print "chr"$2"_"$3"_"$4"_"$5, $0 }' tmp.txt \
				> with_unique_id_${PHENO}_${anc}_${ss}.txt
				rm tmp.txt

				# Create bed file with chrX, pos0, pos1, unique_id
				awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$3, $4-1, $4, $1}' \
				with_unique_id_${PHENO}_${anc}_${ss}.txt \
				> ${PHENO}_${anc}_${ss}_hg19.bed
					
				liftOver ${PHENO}_${anc}_${ss}_hg19.bed \
				/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
				${PHENO}_${anc}_${ss}_hg38.bed \
				${PHENO}_${anc}_${ss}_unlifted.bed

				# Log the % of SNPs succesfully lifted over
				hg19=$(cat ${PHENO}_${anc}_${ss}_hg19.bed | wc -l)
				printf "\t # SNPs in hg19: $hg19 \n" >> lifted/liftover_hg19Tohg38.log
				hg38=$(cat ${PHENO}_${anc}_${ss}_hg38.bed | wc -l)
				printf "\t # SNPs lifted over to hg38: $hg38 \n" >> lifted/liftover_hg19Tohg38.log
				perc_lifted=$(echo "scale=3;${hg38}/${hg19}" | bc)
				printf "\t percent successfully lifted: $(echo "scale=3;${hg38}/${hg19}*100" | bc) \n" >> lifted/liftover_hg19Tohg38.log

				# Create a new sumstats file by joining via unique_id
				# unique_id: 4th column in the lifted over bed file and 1st column in the original sumstats
				join -1 4 -2 1 -o 2.1,1.1,1.3,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 \
				<(sort -k 4 ${PHENO}_${anc}_${ss}_hg38.bed) <(sort -k 1 with_unique_id_${PHENO}_${anc}_${ss}.txt) \
				> lifted/tmp.txt				
				# Add in ID by pasting chr to hg38 pos
				awk -v FS=" " -v OFS=" " '{ print $1,$2":"$3, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' lifted/tmp.txt \
				> lifted/${PHENO}_${anc}_${ss}_hg38.txt
				sed -i '1s/^/chrCHROM_GENPOS_ALLELE1_ALLELE0 ID CHROM GENPOS ALLELE1 ALLELE0 A1FREQ MAF BETA SE CHISQ N PVALUE\n/' \
				lifted/${PHENO}_${anc}_${ss}_hg38.txt
				rm lifted/tmp.txt
				mv *.bed lifted/tmp_files
			fi
		done
	done
done 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

