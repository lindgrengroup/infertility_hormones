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
grep -w 'K5B1.\|K5B1z\|C139.\|X40L4\|C1341\|K5B1z\|K5B10\|C1341\|C1342\|X401I\|C1322\|C1323\|XE2tA\|C1321\|C133.\|XE2Q2\|C1330\|X40L0\|C133y\|XE10V\|C133z\|C1340\|C1343\|C1343\|C1344\|XE10W\|C1320\|C1320\|C132z\|C132.\|X40Ku\|C132y\|C132.\|C1344\|XE10S\|XE10T\|X40LI\|X40L9\|Xa9Ap' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > ukbb_participants_with_IHH_readcodes.txt

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > header_gp_clinical.txt

cat header_gp_clinical.txt ukbb_participants_with_IHH_readcodes.txt > ukbb_participants_with_IHH_readcodes.tmp

mv ukbb_participants_with_IHH_readcodes.tmp ukbb_participants_with_IHH_readcodes.txt

rm header_gp_clinical.txt


