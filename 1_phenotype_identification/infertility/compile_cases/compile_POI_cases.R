setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")


# ICD-code based cases: 82 unique cases (pre-exclusions)
icd<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_POI_cases_icd9_10.txt",sep="\t",header=T,stringsAsFactors=F)
icd_POI1<-icd[which(icd$diag_icd10 %in% c("E283","E288","E289")),]
icd_POI2<-icd[which(icd$diag_icd9 %in% c("2563","2568","2569")),]
icd_POI<-rbind(icd_POI1,icd_POI2)
icd_POI<-icd_POI[,which(names(icd_POI) %in% c("eid","diag_icd9","diag_icd10"))]
icd_POI<-icd_POI[duplicated(icd_POI$eid)==F,]
names(icd_POI)[c(2,3)]<-c("diag_icd9_POI","diag_icd10_POI")

icd10_exclusions<-icd[which(icd$diag_icd10 %in% 
c("Q96","E894","E282","E231","E893","E24","D352","D353","C751","C752","E22","N643","D44","C56","C74","D391","Q564","B20","B21","B22","B23","B24","Z21","E230","E233","E236","E237","R630")),]
icd10_exclusions<-icd10_exclusions[,which(names(icd10_exclusions) %in% c("eid","diag_icd10"))]
icd10_exclusions<-icd10_exclusions[duplicated(icd10_exclusions$eid)==F,]
names(icd10_exclusions)[2]<-"diag_icd10_excl"

icd9_exclusions<-icd[which(icd$diag_icd9 %in% 
c("7586","2562","2564","2537","2550","2273","1943","2530","2531","2536","6116","2374","2372","2370","2371","2373","2580","1830","1940","2362","7527","42","V08","2532","2538","2539","7830")),]
icd9_exclusions<-icd9_exclusions[,which(names(icd9_exclusions) %in% c("eid","diag_icd9"))]
icd9_exclusions<-icd9_exclusions[duplicated(icd9_exclusions$eid)==F,]
names(icd9_exclusions)[2]<-"diag_icd9_excl"

icd<-merge(icd_POI,icd10_exclusions,by="eid",all.x=T,all.y=T)
icd<-merge(icd,icd9_exclusions,by="eid",all=T)

# Primary care based cases.  207 cases based on read2 and 391 based on read3
pc<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/primary_care/ukbb_participants_with_POI_readcodes.txt",sep="\t",header=T,stringsAsFactors=F)
pc<-pc[,which(names(pc) %in% c("eid","read_2","read_3"))]

read2_rows<-which(pc$read_2 %in% c("C1633","C1630","C1631","C1634","C163z","C163y","C1632","C163."))
read3_rows<-which(pc$read_3 %in% c("C1630","XE10k","X408w","C163z","C163y","C1632","C163.","X40NH","X40NK","XSESu","X40NJ","X40NG","X408v","XE10j","X408W","XaZ6K"))

read2_case<-pc[read2_rows,]
read2_case<-read2_case[,-which(names(read2_case) %in% "read_3")]
names(read2_case)<-c("eid","read2_code_case")
read2_case<-read2_case[duplicated(read2_case$eid)==F,]

read3_case<-pc[read3_rows,]
read3_case<-read3_case[,-which(names(read3_case) %in% "read_2")]
names(read3_case)<-c("eid","read3_code_case")
read3_case<-read3_case[duplicated(read3_case$eid)==F,]


read2_excl_rows<-which(pc$read_2 %in% c("B5420","B542z","B5421","B7H22","B7H2z","B7H2.","B7H20","B7H21","B7H23","B912.","C1320","C132y","C132.","C1323","C1341","C132z","C134.","C1330","C1344","K5B1z","C133.","K5B10","K5B11","C1343","C1321","C1322","C139.","C134z","C133y","C133z","C1340","C1342","K5B1.","C1370","C137.","Cyu4N","C13X.","C136z","C13yz","C13y2","C136."
,"C13y0","C13A.","C13y1","C13y.","C13y4","Cyu44","C13z.","C152.","C164.","C165.","C137z","C1371","C1372","C1620","C162z","C1621","C162.","C1622","C162y","K316.","PC7..","PC7z0","PC7z.","R030z","R030."))
read3_excl_rows<-which(pc$read_3 %in% c("B542z","B7H22","B7H2z","XE2vb","B7H20","X78aT","B912.","C1320","C132y","XE10Q","C1341","C132z","C134."
,"C1330","C1344","K5B1z","XE2Q2","K5B10","K5B11","C1343","C1321","C1322","XE10X","XE10V","C133z","C1340","C1342"
,"K5B1.","C1370","C137.","Cyu4N","C136z","C13yz","C13y2","C136.","C13y0","C13y1","C13y.","C13y4","Cyu44",	
"C13z.","XE10b","XE10l","X406n","C137z","C1371","C1372","C1620","C162z","C1621","C162.","C1622","C162y",	
"K316.","PC7..","PC7z0","PC7z.","R030z","R030.","X70M6","A7881","A7885","X20Il","XaFrh","X78aB","B5420",	
"B5421","X77nb","X78aO","XC0u9","X78aP","X77nd","X77nc","X78aI","X78aQ","X78aE","X78aL","X78aM","X78aN",	
"X78aJ","X78aK","X78aD","B7H2.","X406p","X40NU","X78Wn","X406o","X78X5","X78Wm","X78X1","X78dr","B9...",
"XE10N","X40L0","XE10T","XE10S","X408Y","X40Ky","X40L6","XE2tA","X40Or","X40L1","X40LG","X40Ku","X40L8",
"X401I","XE10R","X40Kt","X40Kx","Xa8A9","X40Ks","C133y","X40LI","X40L3","X40L4","X40Km","Xa9Ap","C134z",	
"XE10W","X40L2","X40LE","X40L7","X40Kw","X40Kn","C132.","X40Kv","X00FL","X40L9","XE10Z","X40LF","X40LH",	"X40Lk",	"X40Ls",	"X40LL",	"X40LP",	"X40LQ",	"X40M9",	"X40Lc",	"X40LO",	"X40MA",	"X40Lj",	"X40LM",	"X40LV",	"X40Ld",	"X40Lb",	"X40LW",	"X40La",	"X40Lt",	"X40LN","X40M4",	"X40LX",	"X40Lu","C150.","X40Kf","X40NW","X40NO","X406j","C164.","X40NE","XaZPs","X40Nd","XM0B1","XE0ew","X785Q",
"XE1LI","Xa22s","XE1LJ","PJ...","X78Ex","XE24f","X76cG","14N8.","XE0Gn","ZV6G7","XE0ea","K5333"))

read2_excl<-pc[read2_excl_rows,]
read2_excl<-read2_excl[,-which(names(read2_excl) %in% "read_3")]
names(read2_excl)<-c("eid","read2_code_excl")
read2_excl<-read2_excl[duplicated(read2_excl$eid)==F,]

read3_excl<-pc[read3_excl_rows,]
read3_excl<-read3_excl[,-which(names(read3_excl) %in% "read_2")]
names(read3_excl)<-c("eid","read3_code_excl")
read3_excl<-read3_excl[duplicated(read3_excl$eid)==F,]

pc<-merge(read2_case,read3_case,by="eid",all=T)
pc<-merge(pc,read2_excl,by="eid",all=T)
pc<-merge(pc,read3_excl,by="eid",all=T)
pc$pc_case<-ifelse(!is.na(pc$read2_code_case) | !is.na(pc$read3_code_case),1,0)
pc$pc_exclusion<-ifelse(!is.na(pc$read2_code_excl) | !is.na(pc$read3_code_excl),1,0)

# Merge in all UKBB participant IDs, age and sex. 
ukbb<-read.table("/well/lindgren/UKBIOBANK/laura/bfat_ffmi/anthro_traits_ukbb.txt",sep="\t",header=T,stringsAsFactors=F)
ukbb<-ukbb[,which(names(ukbb) %in% c("eid","sex","age_at_recruitment"))]

ukbb<-merge(ukbb,icd,by="eid",all.x=T,all.y=T)
ukbb<-merge(ukbb,pc,by="eid",all.x=T,all.y=T)

# Drop participants who have withdrawn consent. 209 participants removed. 
drop<-read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/w11867_20220525.csv",sep=",",header=F,stringsAsFactors=F)
names(drop)[1]<-"eid"
ukbb<-ukbb[-which(ukbb$eid %in% drop$eid),]

# Only keep women
ukbb<-ukbb[ukbb$sex==0,]

# Identify POI cases. 614 POI cases without diagnosis on exclusion criteria list. (670 with exclusion diagnosis)
ukbb$POI_case<-ifelse(!is.na(ukbb$diag_icd10_POI) | !is.na(ukbb$diag_icd9_POI) | (!is.na(ukbb$pc_case) & ukbb$pc_case==1),1,0)
ukbb$exclusion<-ifelse(!is.na(ukbb$diag_icd10_excl) | !is.na(ukbb$diag_icd9_excl) | (!is.na(ukbb$pc_exclusion) & ukbb$pc_exclusion==1),1,0)
ukbb$POI_case_include<-ifelse(ukbb$POI_case==1 & ukbb$exclusion==0,1,0)

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(ukbb,"ukbb_participants_POI_cases_exclusions.txt",sep="\t",quote=F,row.names=F)

#

input<-ukbb[,c(1,1,16)]
names(input)<-c("FID","IID","POI")
write.table(input,"/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/POI_eur_regenie.txt",sep="\t",quote=F,row.names=F)

# Drop non-Europeans from pheno input file
data<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/POI_eur_regenie.txt",sep="\t",header=T,stringsAsFactors=F)

eur<-read.table("/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt",header=T,sep="\t",stringsAsFactors=F)
eur<-eur[eur$genetically_european==1,]

data<-data[which(data$IID %in% eur$eid),]

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(data,"POI_eur_regenie_europeans_only.txt",sep="\t",quote=F,row.names=F)



