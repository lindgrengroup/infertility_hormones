setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")


# ICD-code based cases: 532 unique cases (pre-exclusions)
icd<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_IHH_cases_icd9_10.txt",sep="\t",header=T,stringsAsFactors=F)
icd_IHH<-icd[icd$diag_icd10=="E230",]
icd_IHH<-icd_IHH[,which(names(icd_IHH) %in% c("eid","diag_icd10"))]
icd_IHH<-icd_IHH[duplicated(icd_IHH$eid)==F,]
names(icd_IHH)[2]<-"diag_icd10_IHH"
icd_exclusions<-icd[which(icd$diag_icd10 %in% 
c("E231","E232","E893","E221","E8311","E24","D0Y07ZZ","D0Y17ZZ","D0Y0FZZ","D0Y1FZZ","D352",
"D444","C752","D353","E220","E228","E229")),]
icd_exclusions<-icd_exclusions[,which(names(icd_exclusions) %in% c("eid","diag_icd10"))]
icd_exclusions<-icd_exclusions[duplicated(icd_exclusions$eid)==F,]
names(icd_exclusions)[2]<-"diag_icd10_excl"

icd<-merge(icd_IHH,icd_exclusions,by="eid",all.x=T,all.y=T)


# Primary care based cases.  37 cases based on read2 and 69 based on read3
pc<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/primary_care/ukbb_participants_with_IHH_readcodes.txt",sep="\t",header=T,stringsAsFactors=F)
pc<-pc[,which(names(pc) %in% c("eid","read_2","read_3"))]

read2_rows<-which(pc$read_2 %in% c("K5B1.","C139.","K5B1z","K5B10","C1341","C1342"))
read3_rows<-which(pc$read_3 %in% c("K5B1z","X40L4","C1341","X401I"))

read2_case<-pc[read2_rows,]
read2_case<-read2_case[,-which(names(read2_case) %in% "read_3")]
names(read2_case)<-c("eid","read2_code_case")
read2_case<-read2_case[duplicated(read2_case$eid)==F,]

read3_case<-pc[read3_rows,]
read3_case<-read3_case[,-which(names(read3_case) %in% "read_2")]
names(read3_case)<-c("eid","read3_code_case")
read3_case<-read3_case[duplicated(read3_case$eid)==F,]


read2_excl_rows<-which(pc$read_2 %in% c("C1322","C1323","C1321","C133.","C1330","C133y","C133z","C1340","C1343", "C1344","C1320","C132z","C132.","C132y"))
read3_excl_rows<-which(pc$read_3 %in% c("XE2tA","XE2Q2","X40L0","XE10V","C1343","XE10W","C1320","X40Ku","C132.","C1344",
"XE10S","XE10T","X40LI","X40L9","Xa9Ap"))

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

# Identify IHH cases. 303 IHH cases without diagnosis on exclusion criteria list. (611 with exclusion diagnosis)
ukbb$IHH_case<-ifelse(!is.na(ukbb$diag_icd10_IHH) | (!is.na(ukbb$pc_case) & ukbb$pc_case==1),1,0)
ukbb$exclusion<-ifelse(!is.na(ukbb$diag_icd10_excl) | (!is.na(ukbb$pc_exclusion) & ukbb$pc_exclusion==1),1,0)
ukbb$IHH_case_include<-ifelse(ukbb$IHH_case==1 & ukbb$exclusion==0,1,0)

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(ukbb,"ukbb_participants_IHH_cases_exclusions.txt",sep="\t",quote=F,row.names=F)

#

input<-ukbb[,c(1,1,14)]
names(input)<-c("FID","IID","IHH")
write.table(input,"/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/ihh_eur_regenie.txt",sep="\t",quote=F,row.names=F)


# Drop non-Europeans from pheno input file
data<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/ihh_eur_regenie.txt",sep="\t",header=T,stringsAsFactors=F)

eur<-read.table("/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt",header=T,sep="\t",stringsAsFactors=F)
eur<-eur[eur$genetically_european==1,]

data<-data[which(data$IID %in% eur$eid),]

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(data,"ihh_eur_regenie_europeans_only.txt",sep="\t",quote=F,row.names=F)


