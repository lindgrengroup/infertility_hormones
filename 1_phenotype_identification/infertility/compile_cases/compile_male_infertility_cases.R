setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")

# Self-reported male infertility cases. 36 unique cases. 
selfrep<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/self_reported/male_infertility_cases_self_rep_ukbb.txt", header=F,stringsAsFactors=F)
selfrep<-selfrep[,-3]
selfrep$self_reported<-1
names(selfrep)[1:2]<-c("eid","self_rep_code")
selfrep<-selfrep[duplicated(selfrep$eid)==F,]

# ICD-code based cases
prim<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_male_infertility_cases_icd9_10_primary.txt",sep="\t",header=T,stringsAsFactors=F)
sec10<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_male_infertility_cases_icd10_secondary.txt",sep="\t",header=T,stringsAsFactors=F)
sec9<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_male_infertility_cases_icd9_secondary.txt",sep="\t",header=T,stringsAsFactors=F)

# Primary ICD diagnoses - none for ICD-9, 63 unique cases based on ICD10
prim<-prim[,which(names(prim) %in% c("eid","diag_icd10","diag_icd9"))]
prim<-prim[which(prim$diag_icd10 %in% c("N46")),]
prim$prim_ICD10<-1
prim<-prim[,-which(names(prim) %in% "diag_icd9")]
names(prim)<-c("eid","prim_ICD10_code","prim_ICD10")
prim<-prim[duplicated(prim$eid)==F,]

# Secondary ICD-10 diagnoses - 13 unique cases
sec10<-sec10[which(sec10$diag_icd10 %in% c("N46")),]
sec10<-sec10[,which(names(sec10) %in% c("eid","diag_icd10"))]
sec10$sec_ICD10<-1
names(sec10)<-c("eid","sec_ICD10_code","sec_ICD10")
sec10<-sec10[duplicated(sec10$eid)==F,]

# Primary care based cases. 1403 read3 cases, 543 read2 cases. 
pc<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/primary_care/ukbb_participants_with_male_infertility_readcodes.txt",sep="\t",header=T,stringsAsFactors=F)
pc<-pc[,which(names(pc) %in% c("eid","read_2","read_3"))]

read2_rows<-which(pc$read_2 %in% c("4914","4915","4916","4917","K26..","K260.","K261.","K2610","K262.","K26y.","K26y2","K26yz","K26z.","K5B6.","K26y4","K5B10","K5B11","K5By0","K5By1","K5Byz"))
read3_rows<-which(pc$read_3 %in% c("4914","4915","4916","4917","K26..","XE0eA","K261.","XaXbT","XaXUx","K26y.","K26yz","K26z.","X409B","K26y4","K5B10","K5B11","K5By0","K5By1","K5Byz"))

read2<-pc[read2_rows,]
read2$read2<-1
read2<-read2[,-which(names(read2) %in% "read_3")]
names(read2)<-c("eid","read2_code","read2")
read2<-read2[duplicated(read2$eid)==F,]

read3<-pc[read3_rows,]
read3$read3<-1
read3<-read3[,-which(names(read3) %in% "read_2")]
names(read3)<-c("eid","read3_code","read3")
read3<-read3[duplicated(read3$eid)==F,]


# Merge in all UKBB participant IDs, age and sex. 
ukbb<-read.table("/well/lindgren/UKBIOBANK/laura/bfat_ffmi/anthro_traits_ukbb.txt",sep="\t",header=T,stringsAsFactors=F)
ukbb<-ukbb[,which(names(ukbb) %in% c("eid","sex","age_at_recruitment"))]

ukbb<-merge(ukbb,selfrep,by="eid",all.x=T,all.y=T)
ukbb<-merge(ukbb,prim,by="eid",all.x=T,all.y=T)
ukbb<-merge(ukbb,sec10,by="eid",all.x=T,all.y=T)
ukbb<-merge(ukbb,read3,by="eid",all.x=T,all.y=T)
ukbb<-merge(ukbb,read2,by="eid",all.x=T,all.y=T)

# Drop participants who have withdrawn consent. 209 participants removed. 
drop<-read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/w11867_20220525.csv",sep=",",header=F,stringsAsFactors=F)
names(drop)[1]<-"eid"
ukbb<-ukbb[-which(ukbb$eid %in% drop$eid),]

# Drop women
ukbb<-ukbb[ukbb$sex==1,]

# Identify male infertility cases. 666 female infertility cases
ukbb$male_infertility<-ifelse(!is.na(ukbb$self_reported) | !is.na(ukbb$prim_ICD10) | !is.na(ukbb$sec_ICD10) | !is.na(ukbb$read3) | !is.na(ukbb$read2),1,0)

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(ukbb,"male_ukbb_participants_infertility_cases.txt",sep="\t",quote=F,row.names=F)

#

input<-ukbb[,c(1,1,14)]
names(input)<-c("FID","IID","male_infertility")
write.table(input,"/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/male_infertility_traits_eur_regenie.txt",sep="\t",quote=F,row.names=F)

# Drop non-Europeans from pheno input file
data<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/male_infertility_traits_eur_regenie.txt",sep="\t",header=T,stringsAsFactors=F)

eur<-read.table("/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt",header=T,sep="\t",stringsAsFactors=F)
eur<-eur[eur$genetically_european==1,]

data<-data[which(data$IID %in% eur$eid),]

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(data,"male_infertility_traits_eur_regenie_europeans_only.txt",sep="\t",quote=F,row.names=F)



