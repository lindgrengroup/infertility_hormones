#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 31/10/22

# LiftOver public sumstats from hg19 to hg38

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J liftover_public_hg19_to_hg38
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/liftover_public_hg19_to_hg38-%j.out

echo "########################################################"
echo "Slurm Job ID: ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load liftOver/20210519

for PHENO in FSH LH Oestradiol Progesterone Testosterone; do
	cd /well/lindgren/samvida/hormones_infertility/public_sumstats/${PHENO}
	mkdir lifted
	mkdir lifted/tmp_files
	touch lifted/liftover_hg19Tohg38.log

	cat strata_list.txt | while read STRATA; do
		printf "** Strata: ${STRATA}\n" >> lifted/liftover_hg19Tohg38.log

		awk -v FS="\t" -v OFS="\t" '{ sub("23", "X", $2) } 1' filtered_hg19/${STRATA}_filtered.txt > tmp.txt
		# Add unique_identifier to sumstats (chr_pos_A1_A2)
		awk -v FS="\t" -v OFS="\t" '{ print "chr"$2"_"$3"_"$4"_"$5, $0 }' tmp.txt \
		> with_unique_id_${STRATA}_filtered.txt
		rm tmp.txt
		
		# Create bed file with chrX, pos0, pos1, unique_id (CHR_POS_A1_A2)
		awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$3, $4-1, $4, $1}' \
		with_unique_id_${STRATA}_filtered.txt \
		> ${STRATA}_hg19.bed
				
		liftOver ${STRATA}_hg19.bed \
		/well/lindgren/samvida/Resources/hg19ToHg38.over.chain.gz \
		${STRATA}_hg38.bed \
		${STRATA}_unlifted.bed

		# Log the % of SNPs succesfully lifted over
		hg19=$(cat ${STRATA}_hg19.bed | wc -l)
		printf "\t # SNPs in hg19: $hg19 \n" >> lifted/liftover_hg19Tohg38.log
		hg38=$(cat ${STRATA}_hg38.bed | wc -l)
		printf "\t # SNPs lifted over to hg38: $hg38 \n" >> lifted/liftover_hg19Tohg38.log
		perc_lifted=$(echo "scale=3;${hg38}/${hg19}" | bc)
		printf "\t percent successfully lifted: $(echo "scale=3;${hg38}/${hg19}*100" | bc) \n" >> lifted/liftover_hg19Tohg38.log

		# Create a new sumstats file by joining via the unique_id
		# unique_id: 4th column in the lifted over bed file and 1st column in the original sumstats
		# Add header
		join -1 4 -2 1 -o 2.1,1.1,1.3,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 \
		<(sort -k 4 ${STRATA}_hg38.bed) <(sort -k 1 with_unique_id_${STRATA}_filtered.txt) \
		> lifted/tmp.txt
		# Add in ID by pasting chr to hg38 pos
		awk -v FS=" " -v OFS=" " '{ print $1, $2":"$3, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' lifted/tmp.txt \
		> lifted/${STRATA}_hg38.txt
		sed -i '1s/^/chrCHROM_GENPOS_ALLELE1_ALLELE0 ID CHROM GENPOS ALLELE1 ALLELE0 A1FREQ MAF BETA SE STAT N PVALUE\n/' \
		lifted/${STRATA}_hg38.txt
		mv *.bed lifted/tmp_files
	done 
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
