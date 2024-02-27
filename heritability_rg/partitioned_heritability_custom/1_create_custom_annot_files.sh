#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 16/03/23

# LD-scores from LDSC Bulik-Sullivan et al.

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J create_custom_annot_files
#SBATCH -o /well/lindgren/samvida/hormones_infertility/ovary_celltype_markers_melody/logs/create_custom_annot_files-%j.out

echo "########################################################"
echo "Slurm Job ID: ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

# # Create gene-coordinates file
# awk -F'\t' -v OFS='\t' '{print $4, $1, $2, $3}' /well/lindgren/samvida/Resources/hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed \
# > /well/lindgren/samvida/Resources/hg38/hgnc_gene_coordinates.txt
# sed -i '1s/^/GENE\tCHR\tSTART\tEND\n/' /well/lindgren/samvida/Resources/hg38/hgnc_gene_coordinates.txt

# # Create cluster-gene files
# cd /well/lindgren/samvida/hormones_infertility/ovary_celltype_markers_melody

# while IFS=$'\t' read -r -a CLUSTER_GENESET
# do
# 	awk 'NR > 1 {print $1}' ${CLUSTER_GENESET[0]} > for_ldsc/${CLUSTER_GENESET[1]}_hgnc.txt
# done < all_filenames.txt 

cd /well/lindgren/samvida/hormones_infertility/ovary_celltype_markers_melody

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

source /well/lindgren/samvida/python/ldsc-${MODULE_CPU_TYPE}/bin/activate
module load BEDTools/2.30.0-GCC-10.3.0

for CHR in {1..22}; do
	# Make annot files
	python /well/lindgren/samvida/Resources/LDSC/custom_make_annot.py \
	--gene-set-file for_ldsc/${CLUSTER_GENESET}_hgnc.txt \
	--gene-coord-file /well/lindgren/samvida/Resources/hg38/hgnc_gene_coordinates.txt \
	--windowsize 100000 \
	--bimfile /well/lindgren/samvida/Resources/1000Genomes/chr_bedfiles_with_rsid_cm/1000G_Phase3_chr${CHR}.bim \
	--annot-file for_ldsc/annot_files/${CLUSTER_GENESET}.chr${CHR}.annot.gz
done

deactivate
module purge

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

for CHR in {1..22}; do
	ldsc.py \
	--l2 \
	--bfile /well/lindgren/samvida/Resources/1000Genomes/chr_bedfiles_with_rsid_cm/1000G_Phase3_chr${CHR} \
	--ld-wind-cm 1 \
	--annot for_ldsc/annot_files/${CLUSTER_GENESET}.chr${CHR}.annot.gz \
	--thin-annot \
	--out for_ldsc/annot_files/${CLUSTER_GENESET}.chr${CHR} \
	--print-snps /well/lindgren/samvida/Resources/HAPMAP3/hm3_chr${CHR}_only_snp.snplist
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
