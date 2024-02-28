#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_endo

#SBATCH -o ukbb_endo-%j.out 

#SBATCH -e ukbb_endo-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
 


cd /well/lindgren/UKBIOBANK/laura/infertility/primary_care



# Extract participants with read codes for female infertility. 

grep -w 'BBL1.\|BBL2.\|fr...\|fr4..\|K50..\|K500.\|K5000\|K5002	K500z\|K501.\|K502.\|K503.\|K5030\|K5031\|K5032\|K5033\|K503z\|K504.\|K5040\|K5041\|K504z\|K505.\|K5050\|K5051\|K5052\|K505z\|K506.\|K50y.\|K50y0\|K50y1\|K50y2\|K50y3\|K50yz\|K50z.\|Kyu90\|X101N\|X408N\|X408O\|X408T\|X408V\|XaEWA\|XaKcc\|XC04O\|XE0eW\|XE0eX\|7E0D8\|K5001' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > ukbb_participants_with_endometriosis_readcodes.txt

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > header_gp_clinical.txt

cat header_gp_clinical.txt ukbb_participants_with_endometriosis_readcodes.txt > ukbb_participants_with_endometriosis_readcodes.tmp

mv ukbb_participants_with_endometriosis_readcodes.tmp ukbb_participants_with_endometriosis_readcodes.txt

rm header_gp_clinical.txt


