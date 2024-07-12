#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 12/07/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J prepare_loci_file
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/prep_loci_files-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/samvida/hormones_infertility/lava_local/locus_files

# Use awk to add a header "locus" and a sequential count column
awk 'BEGIN { OFS="\t" } 
     NR==1 { print "LOCUS", "CHR", "START", "STOP" } 
     NR>1 { print NR-1, $0 }' LAVA_s2500_m25_f1_w200.blocks > LAVA_s2500_m25_f1_w200.loci

# Split into a separate file per chromosome
# Save the header
header=$(head -n 1 LAVA_s2500_m25_f1_w200.loci)

# Loop through the chromosomes
for chrom in {1..22}; do
    # Print header to the output file
    echo -e "$header" > LAVA_s2500_m25_f1_w200_chr${chrom}.loci
    # Use awk to print lines where the CHR column matches the loop index
    awk -v chrcheck="$chrom" '$2 == chrcheck' LAVA_s2500_m25_f1_w200.loci >> LAVA_s2500_m25_f1_w200_chr${chrom}.loci
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
