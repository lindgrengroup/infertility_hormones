#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_uf

#SBATCH -o ukbb_uf-%j.out 

#SBATCH -e ukbb_uf-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 

# Extract self-reported uterine leiomyoma cases from UKBB (1351)

cd /well/lindgren/UKBIOBANK/laura/infertility/self_reported

for i in {1522..1608};
	do
	VAR_NAME=$(head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab | awk -v var="$i" '{print $var}')
	#echo $VAR_NAME
	cat /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab | awk -v var="$i" '$var == 1351' | awk -v var="$i" -v var2="$VAR_NAME" '{print $1, $var, var2}' >> uterine_leiomyoma_cases_self_rep_ukbb.txt
done

