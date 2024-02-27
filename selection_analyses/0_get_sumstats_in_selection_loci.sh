#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 02/10/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J get_sumstats_selection_loci
#SBATCH -o /well/lindgren/samvida/hormones_infertility/selection_analyses/logs/get_selection_loci-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren/samvida/hormones_infertility/selection_analyses

# Sort hormone and infertility sumstats to start with
# For hormones
while IFS=$'\t' read -r -a STRATA_INFO
do
	# Hormone file columns: ID, CHROM, GENPOS, MAF, Allele1, Allele2, Freq1, FreqSE, BETA, SE, PVALUE, Direction, HetPVal
	# Sort the input 
	sort -k2,2n -k3,3n ${STRATA_INFO[1]} > data/tmp_sorted_${STRATA_INFO[0]}.txt
done < hormone_files_list.txt
# For infertility
while IFS=$'\t' read -r -a STRATA_INFO
do
	# Infertility file columns: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, P.value, HetPVal, N_CASES, N_CONTROLS, maf, chr, pos
	# Sort the input 
	sort -k12,12n -k13,13n ${STRATA_INFO[1]} > data/tmp_sorted_${STRATA_INFO[0]}.txt
done < infertility_files_list.txt

# Subset summary statistics (not rare, MAF > 1%) in selection loci
while IFS=$'\t' read -r -a LOCUS
do
	# For hormones
	while IFS=$'\t' read -r -a STRATA_INFO
	do
		# Hormone file columns: ID, CHROM, GENPOS, MAF, Allele1, Allele2, Freq1, FreqSE, BETA, SE, PVALUE, Direction, HetPVal
		# Intersect with genomic interval
		awk -v CHR="${LOCUS[2]}" -v POS1="${LOCUS[5]}" -v POS2="${LOCUS[6]}" '{ if ($2 == CHR && $3 >= POS1 && $3 <= POS2 && $4 >= 0.01) { print $0 } }' \
		data/tmp_sorted_${STRATA_INFO[0]}.txt > data/chr${LOCUS[2]}_${LOCUS[5]}_${LOCUS[6]}_${STRATA_INFO[0]}.txt
		# Add header column
		sed -i '1s/^/ID\tCHROM\tGENPOS\tMAF\tAllele1\tAllele2\tFreq1\tFreqSE\tBETA\tSE\tPVALUE\tDirection\tHetPVal\n/' \
		data/chr${LOCUS[2]}_${LOCUS[5]}_${LOCUS[6]}_${STRATA_INFO[0]}.txt
	done < hormone_files_list.txt

	# For infertility
	while IFS=$'\t' read -r -a STRATA_INFO
	do
		# Infertility file columns: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, P.value, HetPVal, N_CASES, N_CONTROLS, maf, chr, pos
		# Intersect with genomic interval
		awk -v CHR="${LOCUS[2]}" -v POS1="${LOCUS[5]}" -v POS2="${LOCUS[6]}" '{ if ($12 == CHR && $13 >= POS1 && $13 <= POS2 && $11 >= 0.01) { print $0 } }' \
		data/tmp_sorted_${STRATA_INFO[0]}.txt > data/chr${LOCUS[2]}_${LOCUS[5]}_${LOCUS[6]}_${STRATA_INFO[0]}.txt
		# Add header column
		sed -i '1s/^/MarkerName\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tP.value\tHetPVal\tN_CASES\tN_CONTROLS\tmaf\tchr\tpos\n/' \
		data/chr${LOCUS[2]}_${LOCUS[5]}_${LOCUS[6]}_${STRATA_INFO[0]}.txt
	done < infertility_files_list.txt
done < selection_loci.txt

rm tmp_sorted_*
rm chrCHR_*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0


