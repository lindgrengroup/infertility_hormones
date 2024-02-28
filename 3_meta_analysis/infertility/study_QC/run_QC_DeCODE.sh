#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J decode_QC

#SBATCH --output /well/lindgren/laura/projects/infertility/easyQC/logs/decode_QC-%j.out 

#SBATCH --error /well/lindgren/laura/projects/infertility/easyQC/logs/decode_QC-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=3



echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################
module load R/3.6.2-foss-2019b


while IFS=$'\t' read -r STUDY TRAIT INPUT_FILE N_CASES N_CONTROLS	
	do
	echo $STUDY $TRAIT $INPUT_FILE $N_CASES $N_CONTROLS

LOGFILE=${STUDY}_${TRAIT}_QC.log
Rscript /well/lindgren/laura/projects/infertility/easyQC/make_easyQC_input_DeCODE.R ${STUDY} ${TRAIT} ${INPUT_FILE} ${N_CASES} ${N_CONTROLS} > /well/lindgren/laura/projects/infertility/easyQC/logs/${LOGFILE}
	done < /well/lindgren/laura/projects/infertility/easyQC/DeCODE_QC_inputfile.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

