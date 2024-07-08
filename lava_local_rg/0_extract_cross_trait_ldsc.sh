#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 08/07/2024

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J extract_cross_trait_ldsc
#SBATCH -o /well/lindgren/samvida/hormones_infertility/logs/extract_cross_trait_ldsc-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

# Between hormones and infertility
cd /well/lindgren/samvida/hormones_infertility/rg_hormones_infertility/results

FILES=($(ls rg_EUR_*.log))
N=$(echo ${#FILES[@]})

for I in ${FILES[@]}; do 
    # subset log files to relevant output
	# this will be the last 5 lines (of which the last three lines are about when the analysis ran, so we can remove these)
    tail -n 5 "$I" | head -n 2 > tmp.rg    
    # add to single data set
    if [[ $I == ${FILES[0]} ]]; then
       	cat tmp.rg > /well/lindgren/samvida/hormones_infertility/lava_local/hormones_infertility.rg # only including the header for the first phenotypes
    else
       	cat tmp.rg | sed '1d' >> /well/lindgren/samvida/hormones_infertility/lava_local/hormones_infertility.rg
    fi
	rm tmp.rg
done

# Between hormones/infertility and other reproductive traits
cd /well/lindgren/samvida/hormones_infertility/rg_repro_traits

FILES=($(ls *_vs_repro_phenos.log))
N=$(echo ${#FILES[@]})

for I in ${FILES[@]}; do 
    # subset log files to relevant output
	# this will be the last 11 lines (of which the last three lines are about when the analysis ran, so we can remove these)
    tail -n 11 "$I" | head -n 8 > tmp.rg    
    # add to single data set
    if [[ $I == ${FILES[0]} ]]; then
       	cat tmp.rg > /well/lindgren/samvida/hormones_infertility/lava_local/repro_traits.rg # only including the header for the first phenotypes
    else
       	cat tmp.rg | sed '1d' >> /well/lindgren/samvida/hormones_infertility/lava_local/repro_traits.rg
    fi
	rm tmp.rg
done


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
