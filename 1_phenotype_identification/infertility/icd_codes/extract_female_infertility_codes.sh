#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_icd_female_infertility

#SBATCH -o ukbb_icd_female_infertility-%.out 

#SBATCH -e ukbb_icd_female_infertility-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
 


# From HES based on ICD-9 and ICD-10: female infertility diagnosis
# Primary diagnoses ICD-9 and 10

for ICD in N97 N970 N971 N972 N973 N974 N978 N979 6280 6283 6284 6288 6289
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin.tsv >> ukbb_female_infertility_cases_icd9_10_primary.txt
done


head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin.tsv > header_hes_main.txt
cat header_hes_main.txt ukbb_female_infertility_cases_icd9_10_primary.txt > ukbb_female_infertility_cases_icd9_10_primary.tmp
mv ukbb_female_infertility_cases_icd9_10_primary.tmp ukbb_female_infertility_cases_icd9_10_primary.txt
rm header_hes_main.txt 

# Secondary diagnoses ICD-9
for ICD in 6280 6283 6284 6288 6289
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag9.tsv >> ukbb_female_infertility_cases_icd9_secondary.txt
done

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag9.tsv > header_hes_secondary_icd9.txt
cat header_hes_secondary_icd9.txt ukbb_female_infertility_cases_icd9_secondary.txt > ukbb_female_infertility_cases_icd9_secondary.tmp
mv ukbb_female_infertility_cases_icd9_secondary.tmp ukbb_female_infertility_cases_icd9_secondary.txt
rm header_hes_secondary_icd9.txt

# Secondary diagnoses ICD-10
for ICD in N97 N970 N971 N972 N973 N974 N978 N979 
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag10.tsv >> ukbb_female_infertility_cases_icd10_secondary.txt
done

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/ukb_hesin_diag10.tsv > header_hes_secondary_icd10.txt
cat header_hes_secondary_icd10.txt ukbb_female_infertility_cases_icd10_secondary.txt > ukbb_female_infertility_cases_icd10_secondary.tmp
mv ukbb_female_infertility_cases_icd10_secondary.tmp ukbb_female_infertility_cases_icd10_secondary.txt
rm header_hes_secondary_icd10.txt


