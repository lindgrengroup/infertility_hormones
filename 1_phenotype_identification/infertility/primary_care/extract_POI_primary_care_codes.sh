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
grep -w 'C1633\|C1631\|C1634\|C1630\|XE10k\|C163z\|C163y\|C1632\|C163.\|X40NH\|X40NK\|XSESu\|X40NJ\|X40NG\|X408v\|XE10j\|X408W\|XaZ6K\|B7H21\|B7H23\|C1323\|C133.\|C139.\|C13X.\|C13A.\|C152.\|C165.\|B542z\|B7H22\|B7H2z\|XE2vb\|B7H20\|X78aT\|B912.\|C1320\|C132y\|XE10Q\|C1341\|C132z\|C134.\|C1330\|C1344\|K5B1z\|XE2Q2\|K5B10\|K5B11\|C1343\|C1321\|C1322\|XE10X\|XE10V\|C133z\|C1340\|C1342\|K5B1.\|C1370\|C137.\|Cyu4N\|C136z\|C13yz\|C13y2\|C136.\|C13y0\|C13y1\|C13y.\|C13y4\|Cyu44\|C13z.\|XE10b\|XE10l\|X406n\|C137z\|C1371\|C1372\|C1620\|C162z\|C1621\|C162.\|C1622\|C162y\|K316.\|PC7..\|PC7z0\|PC7z.\|R030z\|R030.\|X70M6\|A7881\|A7885\|X20Il\|XaFrh\|X78aB\|B5420\|B5421\|X77nb\|X78aO\|XC0u9\|X78aP\|X77nd\|X77nc\|X78aI\|X78aQ\|X78aE\|X78aL\|X78aM\|X78aN\|X78aJ\|X78aK\|X78aD\|B7H2.\|X406p\|X40NU\|X78Wn\|X406o\|X78X5\|X78Wm\|X78X1\|X78dr\|B9...\|XE10N\|X40L0\|XE10T\|XE10S\|X408Y\|X40Ky\|X40L6\|XE2tA\|X40Or\|X40L1\|X40LG\|X40Ku\|X40L8\|X401I\|XE10R\|X40Kt\|X40Kx\|Xa8A9\|X40Ks\|C133y\|X40LI\|X40L3\|X40L4\|X40Km\|Xa9Ap\|C134z\|XE10W\|X40L2\|X40LE\|X40L7\|X40Kw\|X40Kn\|C132.\|X40Kv\|X00FL\|X40L9\|XE10Z\|X40LF\|X40LH\|X40Lk\|X40Lr\|X40Ls\|X40LL\|X40LP\|X40LQ\|X40M9\|X40Lc\|X40LO\|X40MA\|X40Lj\|X40LM\|X40LV\|X40Ld\|X40Lb\|X40LW\|X40La\|X40Lt\|X40LN\|X40M4\|X40LX\|X40Lu\|C150.\|X40Kf\|X40NW\|X40NO\|X406j\|C164.\|X40NE\|XaZPs\|X40Nd\|XM0B1\|XE0ew\|X785Q\|XE1LI\|Xa22s\|XE1LJ\|PJ...\|X78Ex\|XE24f\|X76cG\|14N8.\|XE0Gn\|ZV6G7\|X76Wd\|XaCxM\|XE0ea\|K5333' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > ukbb_participants_with_POI_readcodes.txt

head -1 /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt > header_gp_clinical.txt

cat header_gp_clinical.txt ukbb_participants_with_POI_readcodes.txt > ukbb_participants_with_POI_readcodes.tmp

mv ukbb_participants_with_POI_readcodes.tmp ukbb_participants_with_POI_readcodes.txt

rm header_gp_clinical.txt


