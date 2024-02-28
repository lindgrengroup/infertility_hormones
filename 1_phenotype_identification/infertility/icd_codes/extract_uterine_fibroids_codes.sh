#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_icd_uterine_fibroids

#SBATCH -o ukbb_icd_uterine_fibroids-%.out 

#SBATCH -e ukbb_icd_uterine_fibroids-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
 


# From HES based on ICD-9 and ICD-10: uterine_fibroids diagnosis
# Primary diagnoses ICD-9 and 10

for ICD in D25 D250 D251 D252 D259 218
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin.tsv >> ukbb_uterine_fibroids_cases_icd9_10_primary.txt
done


head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin.tsv > header_hes_main.txt
cat header_hes_main.txt ukbb_uterine_fibroids_cases_icd9_10_primary.txt > ukbb_uterine_fibroids_cases_icd9_10_primary.tmp
mv ukbb_uterine_fibroids_cases_icd9_10_primary.tmp ukbb_uterine_fibroids_cases_icd9_10_primary.txt
rm header_hes_main.txt 

# Secondary diagnoses ICD-9
for ICD in 218
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag9.tsv >> ukbb_uterine_fibroids_cases_icd9_secondary.txt
done

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag9.tsv > header_hes_secondary_icd9.txt
cat header_hes_secondary_icd9.txt ukbb_uterine_fibroids_cases_icd9_secondary.txt > ukbb_uterine_fibroids_cases_icd9_secondary.tmp
mv ukbb_uterine_fibroids_cases_icd9_secondary.tmp ukbb_uterine_fibroids_cases_icd9_secondary.txt
rm header_hes_secondary_icd9.txt

# Secondary diagnoses ICD-10
for ICD in D25 D250 D251 D252 D259 
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag10.tsv >> ukbb_uterine_fibroids_cases_icd10_secondary.txt
done

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag10.tsv > header_hes_secondary_icd10.txt
cat header_hes_secondary_icd10.txt ukbb_uterine_fibroids_cases_icd10_secondary.txt > ukbb_uterine_fibroids_cases_icd10_secondary.tmp
mv ukbb_uterine_fibroids_cases_icd10_secondary.tmp ukbb_uterine_fibroids_cases_icd10_secondary.txt
rm header_hes_secondary_icd10.txt


