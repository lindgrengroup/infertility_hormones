setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")

# Self-reported endometriosis cases. 4195 unique cases. 
selfrep<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/self_reported/endometriosis_cases_self_rep_ukbb.txt", header=F,stringsAsFactors=F)
selfrep<-selfrep[,-3]
selfrep$self_reported<-1
names(selfrep)[1:2]<-c("eid","self_rep_code")
selfrep<-selfrep[duplicated(selfrep$eid)==F,]

# ICD-code based cases
prim<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_endometriosis_cases_icd9_10_primary.txt",sep="\t",header=T,stringsAsFactors=F)
sec10<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_endometriosis_cases_icd10_secondary.txt",sep="\t",header=T,stringsAsFactors=F)
sec9<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_endometriosis_cases_icd9_secondary.txt",sep="\t",header=T,stringsAsFactors=F)

# Primary ICD diagnoses - none for ICD-9, 2155 unique cases based on ICD10
prim<-prim[,which(names(prim) %in% c("eid","diag_icd10","diag_icd9"))]
prim<-prim[which(prim$diag_icd10 %in% c("N80", "N800", "N801", "N802","N803","N804" ,"N805","N806","N808","N809")),]
prim$prim_ICD10<-1
prim<-prim[,-which(names(prim) %in% "diag_icd9")]
names(prim)<-c("eid","prim_ICD10_code","prim_ICD10")
prim<-prim[duplicated(prim$eid)==F,]

# Secondary ICD-10 diagnoses - 3249 unique cases
sec10<-sec10[which(sec10$diag_icd10 %in% c("N80", "N800", "N801", "N802","N803","N804" ,"N805","N806","N808","N809")),]
sec10<-sec10[,which(names(sec10) %in% c("eid","diag_icd10"))]
sec10$sec_ICD10<-1
names(sec10)<-c("eid","sec_ICD10_code","sec_ICD10")
sec10<-sec10[duplicated(sec10$eid)==F,]

# Secondary ICD-9 diagngosis - 140 cases
sec9<-sec9[which(sec9$diag_icd9 %in% c("617","6170","6171","6172","6173","6174","6175","6176","6178","6179")),]
sec9<-sec9[,which(names(sec9) %in% c("eid","diag_icd9"))]
sec9$sec_ICD9<-1
names(sec9)<-c("eid","sec_ICD9_code","sec_ICD9")
sec9<-sec9[duplicated(sec9$eid)==F,]

# Primary care based cases. 1969 read3 cases, 1149 read2 cases. 
pc<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/primary_care/ukbb_participants_with_endometriosis_readcodes.txt",sep="\t",header=T,stringsAsFactors=F)
pc<-pc[,which(names(pc) %in% c("eid","read_2","read_3"))]

read2_rows<-which(pc$read_2 %in% c("7E0D8","BBL1.","BBL2.","K50..","K500.","K5000","K5001","K5002","K500z",
"K501.","K502.","K503.","K5030","K5031","K5032","K5033","K503z","K504.","K5040",
"K5041","K504z","K505.","K5050","K5051","K5052","K505z","K506.","K50y.","K50y0",
"K50y1","K50y2","K50y3","K50yz","K50z.","Kyu90"))

read3_rows<-which(pc$read_3 %in% c("BBL1.","BBL2.","fr...","fr4..","K50..","K50..","K500.","K5000","K5002","K500z","K501.","K501.","K502.","K503.","K5030","K5031","K5032","K5033","K503z","K504.","K5040","K5041","K504z","K505.","K5050","K5051","K5052","K505z","K506.","K50y.","K50y0","K50y1","K50y2","K50y3","K50yz","K50z.","Kyu90","X101N","X101N","X408N","X408O","X408T","X408V","XaEWA","XaKcc","XC04O","XE0eW","XE0eX"))

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
ukbb<-merge(ukbb,sec9,by="eid",all.x=T,all.y=T)
ukbb<-merge(ukbb,read3,by="eid",all.x=T,all.y=T)
ukbb<-merge(ukbb,read2,by="eid",all.x=T,all.y=T)

# Drop participants who have withdrawn consent. 209 participants removed. 
drop<-read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/w11867_20220525.csv",sep=",",header=F,stringsAsFactors=F)
names(drop)[1]<-"eid"
ukbb<-ukbb[-which(ukbb$eid %in% drop$eid),]

# Drop men
ukbb<-ukbb[ukbb$sex==0,]


# Identify female infertility cases. 3,870 female infertility cases
ukbb$endometriosis<-ifelse(!is.na(ukbb$self_reported) | !is.na(ukbb$prim_ICD10) | !is.na(ukbb$sec_ICD10) | !is.na(ukbb$sec_ICD9) | !is.na(ukbb$read3) | !is.na(ukbb$read2),1,0)


setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(ukbb,"female_ukbb_participants_endometriosis_diagnoses.txt",sep="\t",quote=F,row.names=F)


















