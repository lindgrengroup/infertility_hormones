setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")

# Self-reported female infertility cases. 8458 unique cases. 
selfrep<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/self_reported/uterine_leiomyoma_cases_self_rep_ukbb.txt", header=F,stringsAsFactors=F)
selfrep<-selfrep[,-3]
selfrep$self_reported<-1
names(selfrep)[1:2]<-c("eid","self_rep_code")
selfrep<-selfrep[duplicated(selfrep$eid)==F,]

# ICD-code based cases
prim<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_uterine_fibroids_cases_icd9_10_primary.txt",sep="\t",header=T,stringsAsFactors=F)
sec10<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_uterine_fibroids_cases_icd10_secondary.txt",sep="\t",header=T,stringsAsFactors=F)
sec9<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_uterine_fibroids_cases_icd9_secondary.txt",sep="\t",header=T,stringsAsFactors=F)

# Primary ICD diagnoses - none for ICD-9, 8685 unique cases based on ICD10
prim<-prim[,which(names(prim) %in% c("eid","diag_icd10","diag_icd9"))]
prim<-prim[which(prim$diag_icd10 %in% c("D25","D250","D251","D252","D259")),]
prim$prim_ICD10<-1
prim<-prim[,-which(names(prim) %in% "diag_icd9")]
names(prim)<-c("eid","prim_ICD10_code","prim_ICD10")
prim<-prim[duplicated(prim$eid)==F,]

# Secondary ICD-10 diagnoses - 5637 unique cases
sec10<-sec10[which(sec10$diag_icd10 %in% c("D25","D250","D251","D252","D259")),]
sec10<-sec10[,which(names(sec10) %in% c("eid","diag_icd10"))]
sec10$sec_ICD10<-1
names(sec10)<-c("eid","sec_ICD10_code","sec_ICD10")
sec10<-sec10[duplicated(sec10$eid)==F,]

# Secondary ICD-9 diagngosis - no cases


# Primary care based cases. 5111 read3 cases, 71795 read2 cases. 
pc<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/primary_care/ukbb_participants_with_UF_readcodes.txt",sep="\t",header=T,stringsAsFactors=F)
pc<-pc[,which(names(pc) %in% c("eid","read_2","read_3"))]

read2_rows<-which(pc$read_2 %in% c("B78..","B780.","B781.","B782.","B78z.","X78Xb","X78XZ","XaNQc"))
read3_rows<-which(pc$read_3 %in% c("B78..","B780.","B781.","B782.","B78z.","X78Xb","X78XZ","XaNQc"))

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

# Drop men
ukbb<-ukbb[ukbb$sex==0,]


# Identify female infertility cases. 3,870 female infertility cases
ukbb$uterine_fibroids<-ifelse(!is.na(ukbb$self_reported) | !is.na(ukbb$prim_ICD10) | !is.na(ukbb$sec_ICD10) | !is.na(ukbb$read3) | !is.na(ukbb$read2),1,0)


setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(ukbb,"female_ukbb_participants_uterine_fibroids_diagnoses.txt",sep="\t",quote=F,row.names=F)

fi<-read.table("female_ukbb_participants_infertility_diagnoses.txt",sep="\t",header=T,stringsAsFactors=F)
fi<-fi[,c(1,16)]

ukbb<-merge(ukbb,fi,by="eid",all=F)

> nrow(ukbb)
[1] 273323
> nrow(ukbb[ukbb$uterine_fibroids==1,])
[1] 24097
> nrow(ukbb[ukbb$female_infertility==1,])
[1] 3870
> nrow(ukbb[ukbb$uterine_fibroids==1 & ukbb$female_infertility==1,])
[1] 628
> nrow(ukbb[ukbb$uterine_fibroids==1 & ukbb$female_infertility==0,])
[1] 23469
> nrow(ukbb[ukbb$uterine_fibroids==0 & ukbb$female_infertility==1,])
[1] 3242
> 628/24097
[1] 0.02606134
>

















