#!/bin/bash 

 
#SBATCH -A lindgren.prj 
#SBATCH -J ukbb_UF

#SBATCH -o ukbb_UF-%j.out 

#SBATCH -e ukbb_UF-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
 


cd /well/lindgren/UKBIOBANK/laura/infertility/primary_care



# Extract participants with read codes for uterine leiomyoma 

grep -w 'B78..\|B780.\|B781.\|B782.\|B78z.\|X78Xb\|X78XZ\|XaNQc' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > ukbb_participants_with_UF_readcodes.txt

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > header_gp_clinical.txt

cat header_gp_clinical.txt ukbb_participants_with_UF_readcodes.txt > ukbb_participants_with_UF_readcodes.tmp

mv ukbb_participants_with_UF_readcodes.tmp ukbb_participants_with_UF_readcodes.txt

rm header_gp_clinical.txt


