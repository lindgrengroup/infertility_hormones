#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J coloc_hormone_inf_F_EUR

#SBATCH --output /well/lindgren/laura/projects/infertility/colocalisations/sex_hormones/logs/coloc_hormone_inf_F_EUR-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/colocalisations/sex_hormones/logs/coloc_hormone_inf_F_EUR-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=4




echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/laura/projects/infertility/colocalisations/sex_hormones

module load R/3.6.2-foss-2019b


for INFERTILITY in  "female_infertility_analysis3_eur" "female_infertility_analysis1_eur" "female_infertility_analysis2_eur"  "female_infertility_analysis4_eur" "female_infertility_analysis5_eur" ;

do 
	for HORMONE in "FSH_F_EUR" "LH_F_EUR" "Oestradiol_F_EUR" "Progesterone_F_EUR" "Testosterone_F_EUR";
		do
			echo "Run Coloc for $INFERTILITY vs $HORMONE"
			Rscript run_coloc_infertility_vs_hormones.R $INFERTILITY $HORMONE
		done
done
