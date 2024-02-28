#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_pcos

#SBATCH -o ukbb_pcos-%j.out 

#SBATCH -e ukbb_pcos-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 



cd /well/lindgren/UKBIOBANK/laura/infertility/self_reported

# Extract PCOS cases self-reported from UKBB (1350)
for i in {1522..1608};

	do
	VAR_NAME=$(head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab | awk -v var="$i" '{print $var}')
	#echo $VAR_NAME
	cat /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab | awk -v var="$i" '$var == 1350' | awk -v var="$i" -v var2="$VAR_NAME" '{print $1, $var, var2}' >> pcos_cases_self_rep_ukbb.txt
	done

