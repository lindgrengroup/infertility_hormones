#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J coloc_hormone_inf_F_all

#SBATCH --output /well/lindgren/laura/projects/infertility/colocalisations/sex_hormones/logs/coloc_hormone_inf_F_all-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/colocalisations/sex_hormones/logs/coloc_hormone_inf_F_all-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=4




echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/laura/projects/infertility/colocalisations/sex_hormones

module load R/3.6.2-foss-2019b


for INFERTILITY in  "female_infertility_analysis3_all" "female_infertility_analysis1_all" "female_infertility_analysis2_all"  "female_infertility_analysis4_all" "female_infertility_analysis5_all" ;

do 
	for HORMONE in "FSH_F_all_anc" "LH_F_all_anc" "Oestradiol_F_all_anc" "Progesterone_F_all_anc" "Testosterone_F_all_anc";
		do
			echo "Run Coloc for $INFERTILITY vs $HORMONE"
			Rscript run_coloc_infertility_vs_hormones.R $INFERTILITY $HORMONE
		done
done
