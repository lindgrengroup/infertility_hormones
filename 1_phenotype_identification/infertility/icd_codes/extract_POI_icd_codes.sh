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

for ICD in 2563 2568 2569 E283 E288 E289 Q96 E894 E282 E231 E893 E24 D352 D353 C751 C752 E22 N643 D44 C56 C74 D391 Q564 B20 B21 B22 B23 B24 Z21 E230 E233 E236 E237 R630 7586 2562 2564 2537 2537 2550 2273 2273 1943 1943 2530 2531 2536 6116 2374 2372 2370 2371 2373 2580 1830 1940 2362 7527 42 V08 2532 2539 2538 2539 7830
do
	grep -w $ICD /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/hesin_diag_21012022.txt >> ukbb_POI_cases_icd9_10.txt
done


head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/HES/hesin_diag_21012022.txt > header_hes_main.txt
cat header_hes_main.txt ukbb_POI_cases_icd9_10.txt > ukbb_POI_cases_icd9_10.tmp
mv ukbb_POI_cases_icd9_10.tmp ukbb_POI_cases_icd9_10.txt
rm header_hes_main.txt 



