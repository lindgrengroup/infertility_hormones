#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_male_inf

#SBATCH -o ukbb_male_inf-%j.out 

#SBATCH -e ukbb_male_inf-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 

# Extract self-reported male infertility cases (1404) from UKBB
cd /well/lindgren/UKBIOBANK/laura/infertility/self_reported

for i in {1522..1608};

	do
	VAR_NAME=$(head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab | awk -v var="$i" '{print $var}')
	#echo $VAR_NAME
	cat /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab | awk -v var="$i" '$var == 1404' | awk -v var="$i" -v var2="$VAR_NAME" '{print $1, $var, var2}' >> male_infertility_cases_self_rep_ukbb.txt
	done

