#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_icd_IHH

#SBATCH -o ukbb_icd_IHH-%.out 

#SBATCH -e ukbb_icd_IHH-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
 


# From HES based on ICD-10: hypogonadism diagnosis
# Primary diagnoses ICD-10

for ICD in E230 E231 E232 E893 E221 E8311 E24 D0Y07ZZ D0Y17ZZ D0Y0FZZ D0Y1FZZ D352 D444 C752 D353 E220 E228 E229
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/hesin_diag_21012022.txt >> ukbb_IHH_cases_icd9_10.txt
done


head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/hesin_diag_21012022.txt > header_hes_main.txt
cat header_hes_main.txt ukbb_IHH_cases_icd9_10.txt > ukbb_IHH_cases_icd9_10.tmp
mv ukbb_IHH_cases_icd9_10.tmp ukbb_IHH_cases_icd9_10.txt
rm header_hes_main.txt 



