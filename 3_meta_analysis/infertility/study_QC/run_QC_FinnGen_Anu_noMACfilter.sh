#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J FinnGen_Anu_QC

#SBATCH --output /well/lindgren/laura/projects/infertility/easyQC/logs/FinnGen_Anu_QC-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/easyQC/logs/FinnGen_Anu_QC-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=6



echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################
module load R/3.6.2-foss-2019b


while IFS=$'\t' read -r STUDY TRAIT INPUT_FILE N_CASES N_CONTROLS	
	do
	echo $STUDY $TRAIT $INPUT_FILE $N_CASES $N_CONTROLS
#	echo "1. Run QC of GWAS summary stats"
	rm /well/lindgren/laura/projects/infertility/easyQC/logs/${STUDY}_${TRAIT}_QC.log
	Rscript /well/lindgren/laura/projects/infertility/easyQC/make_easyQC_input_FinnGen_r10_Anu_noMACfilter.R ${STUDY} ${TRAIT} ${INPUT_FILE} ${N_CASES} ${N_CONTROLS} > /well/lindgren/laura/projects/infertility/easyQC/logs/${STUDY}_${TRAIT}_QC.log

	echo "2. Generate QQ plots"
	Rscript /well/lindgren/laura/projects/infertility/plots/qq_study_level/QQplots_FinnGen_r10_Anu_noMACfilter.R ${STUDY} ${TRAIT}

	echo "3. Generate MH plots"
	Rscript /well/lindgren/laura/projects/infertility/plots/mh_study_level/make_MH_plots_results_Finngen_r10_Anu_noMACfilter.R ${STUDY} ${TRAIT}
	done < /well/lindgren/laura/projects/infertility/easyQC/FinnGen_r10_Anu_female_infertility_QC_inputfile.txt


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

