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

grep -w '4914\|4915\|4916\|4917\|K26..00\|K260.00\|K261.00\|K261000\|K262.00\|K26y.00\|K26y200\|K26yz00\|K26z.00\|K5B6.00\|K26y400\|K5B1000\|K5B1100\|K5By000\|K5By100\|K5Byz11\|K260.11\|K26y300\|X101d\|K26y3\|K26..\|XE0eA\|K261.\|XaXbT\|XaXUx\|K26y.\|K26yz\|K26z.\|X409B\|K26y4\|K5B10\|K5B11\|K5By0\|K5By1\|K5Byz' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > ukbb_participants_with_male_infertility_readcodes.txt

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > header_gp_clinical.txt

cat header_gp_clinical.txt ukbb_participants_with_male_infertility_readcodes.txt > ukbb_participants_with_male_infertility_readcodes.tmp

mv ukbb_participants_with_male_infertility_readcodes.tmp ukbb_participants_with_male_infertility_readcodes.txt

rm header_gp_clinical.txt


