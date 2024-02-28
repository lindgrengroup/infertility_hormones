#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_prim_care

#SBATCH -o ukbb_prim_care-%j.out 

#SBATCH -e ukbb_prim_care-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
 


cd /well/lindgren/UKBIOBANK/laura/infertility/primary_care



# Extract participants with read codes for female infertility. 

grep -w 'K5B..\|K5B0.\|XE0ex\|K5B00\|K5B01\|K5B0z\|K5B1.\|K5B1z\|K5B2.\|K5B20\|K5B21\|X4099\|K5B2z\|K5B3.\|K5B30\|K5B31\|K5B3z\|K5B4.\|K5B40\|K5B41\|K5B5.\|K5B51\|X409B\|XaZ6K\|K5By.\|XE0ey\|K5Bz.\|Kyu9G\|K26y4\|K5B10\|K5B11\|K5By0\|K5By1\|K5Byz\|K5B6.\|K5B7.\|K26y3\|5B23' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > ukbb_participants_with_female_infertility_readcodes.txt

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > header_gp_clinical.txt

cat header_gp_clinical.txt ukbb_participants_with_female_infertility_readcodes.txt > ukbb_participants_with_female_infertility_readcodes.tmp

mv ukbb_participants_with_female_infertility_readcodes.tmp ukbb_participants_with_female_infertility_readcodes.txt

rm header_gp_clinical.txt


