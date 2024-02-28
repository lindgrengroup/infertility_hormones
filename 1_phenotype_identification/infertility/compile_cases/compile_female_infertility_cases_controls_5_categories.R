setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")

# Self-reported female infertility cases. 618 unique cases. 
selfrep<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/self_reported/female_infertility_cases_self_rep_ukbb.txt", header=F,stringsAsFactors=F)
selfrep<-selfrep[,-3]
selfrep$self_reported<-1
names(selfrep)[1:2]<-c("eid","self_rep_code")
selfrep<-selfrep[duplicated(selfrep$eid)==F,]

# ICD-code based cases
prim<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_female_infertility_cases_icd9_10_primary.txt",sep="\t",header=T,stringsAsFactors=F)
sec10<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_female_infertility_cases_icd10_secondary.txt",sep="\t",header=T,stringsAsFactors=F)
sec9<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/icd_codes/ukbb_female_infertility_cases_icd9_secondary.txt",sep="\t",header=T,stringsAsFactors=F)

# Primary ICD diagnoses - none for ICD-9, 1033 unique cases based on ICD10
prim<-prim[,which(names(prim) %in% c("eid","diag_icd10","diag_icd9"))]
prim<-prim[which(prim$diag_icd10 %in% c("N97","N970","N971","N972","N973","N974","N978","N979")),]
prim$prim_ICD10<-1
prim<-prim[,-which(names(prim) %in% "diag_icd9")]
names(prim)<-c("eid","prim_ICD10_code","prim_ICD10")
prim<-prim[duplicated(prim$eid)==F,]

# Secondary ICD-10 diagnoses - 340 unique cases
sec10<-sec10[which(sec10$diag_icd10 %in% c("N97","N970","N971","N972","N973","N974","N978","N979")),]
sec10<-sec10[,which(names(sec10) %in% c("eid","diag_icd10"))]
sec10$sec_ICD10<-1
names(sec10)<-c("eid","sec_ICD10_code","sec_ICD10")
sec10<-sec10[duplicated(sec10$eid)==F,]

# Secondary ICD-9 diagngosis - 30 unique cases
sec9<-sec9[which(sec9$diag_icd9 %in% c(6280,6283,6284,6288,6289)),]
sec9<-sec9[,which(names(sec9) %in% c("eid","diag_icd9"))]
sec9$sec_ICD9<-1
names(sec9)<-c("eid","sec_ICD9_code","sec_ICD9")
sec9<-sec9[duplicated(sec9$eid)==F,]

# Primary care based cases. 1574 read3 cases, 726 read2 cases. 
pc<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/primary_care/ukbb_participants_with_female_infertility_readcodes.txt",sep="\t",header=T,stringsAsFactors=F)
pc<-pc[,which(names(pc) %in% c("eid","read_2","read_3"))]

read2_rows<-which(pc$read_2 %in% c("K5B..","K5B0.","K5B0.","K5B00","K5B01","K5B0z","'K5B1.","K5B1z","K5B2.","K5B20","K5B21","K5B23","K5B2z","K5B3.","K5B30","K5B31","K5B3z","K5B4.","K5B40","K5B41","K5B5.","K5B51","K5B6.","K5B7.","K5By.","K5Byz","K5Bz.","Kyu9G","K26y4","K5B10","K5B11","K5By0","K5By1","K5Byz"))
read3_rows<-which(pc$read_3 %in% c("K5B..","XE0ex","K5B00","K5B01","K5B0z","K5B1.","K5B1z","K5B2.","K5B20","K5B21","X4099","K5B2z","K5B3.","K5B30","K5B31","K5B3z","K5B4.","K5B40","K5B41","K5B5.","K5B51","X409B","XaZ6K","K5By.","XE0ey","K5Bz.","Kyu9G","K26y4","K5B10","K5B11","K5By0","K5By1","K5Byz","K26y3"))

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
ukbb<-merge(ukbb,sec9,by="eid",all.x=T,all.y=T)
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
ukbb$female_infertility<-ifelse(!is.na(ukbb$self_reported) | !is.na(ukbb$prim_ICD10) | !is.na(ukbb$sec_ICD9) | !is.na(ukbb$sec_ICD10) | !is.na(ukbb$read3) | !is.na(ukbb$read2),1,0)

# 499 of 3870 cases only identified through self-report.
nrow(cases[!is.na(cases$self_reported) & is.na(cases$prim_ICD10) & is.na(cases$sec_ICD9) & is.na(cases$sec_ICD10) & is.na(cases$read3) & is.na(cases$read2),])
[1] 499
nrow(cases[!is.na(cases$self_reported),])
[1] 618

# 1978 of 3870 cases only identified through primary care codes. 
nrow(cases[is.na(cases$self_reported) & is.na(cases$prim_ICD10) & is.na(cases$sec_ICD9) & is.na(cases$sec_ICD10) & (!is.na(cases$read3) | !is.na(cases$read2)),])
[1] 1978

# 1146 of 3870 cases only identified through ICD codes. 
nrow(cases[is.na(cases$self_reported) & (!is.na(cases$prim_ICD10) | !is.na(cases$sec_ICD9) | !is.na(cases$sec_ICD10)) & is.na(cases$read3) & is.na(cases$read3),])
[1] 1146

# Breakdown of through which data source infertility cases are identified. 
cases<-ukbb[ukbb$female_infertility==1,]
cases$ICD10<-ifelse(!is.na(cases$prim_ICD10_code),cases$prim_ICD10_code,ifelse(!is.na(cases$sec_ICD10_code),cases$sec_ICD10_code,NA))
nrow(cases[!is.na(cases$ICD10),]) # N=1,297
nrow(cases[!is.na(cases$ICD10) & cases$ICD10=="N970",]) # N=12
nrow(cases[!is.na(cases$ICD10) & cases$ICD10=="N971",]) # N=135
nrow(cases[!is.na(cases$ICD10) & cases$ICD10=="N972",]) # N=8
nrow(cases[!is.na(cases$ICD10) & cases$ICD10=="N973",]) # N=1
nrow(cases[!is.na(cases$ICD10) & cases$ICD10=="N974",]) # N=11
nrow(cases[!is.na(cases$ICD10) & cases$ICD10=="N978",]) # N=20
nrow(cases[!is.na(cases$ICD10) & cases$ICD10=="N979",]) # N=1,110

cases$ICD9<-ifelse(!is.na(cases$prim_ICD10_code),cases$prim_ICD10_code,ifelse(!is.na(cases$sec_ICD10_code),cases$sec_ICD10_code,NA))
nrow(cases[!is.na(cases$ICD10),]) # N=1,297

# Identify anovulatory infertility cases based on ICD10 (N970), ICD9 (6280) and primary care codes (read2 K5B0., K5B00, K5B01, K5B0z; read3: XE0ex, K5B00, K5B01, K5B0z). 54 cases. Mostly from primary care codes (read2 K5B2., K5B20, K5B21, K5B23, K5B2z, K5B3., K5B30, K5B31, K5B3z, K5B4., K5B40, K5B41, K5B5., K5B51). 
ukbb$anovulatory_infertility<-ifelse((!is.na(ukbb$prim_ICD10) & ukbb$prim_ICD10_code=="N970") | (!is.na(ukbb$sec_ICD10) & ukbb$sec_ICD10_code=="N970") | (!is.na(ukbb$sec_ICD9) & ukbb$sec_ICD9_code==6280) | (!is.na(ukbb$read2) & (ukbb$read2_code=="K5B0." | ukbb$read2_code=="K5B00" | ukbb$read2_code=="K5B01" | ukbb$read2_code=="K5B0z")) | (!is.na(ukbb$read3) & (ukbb$read3_code=="XE0ex" | ukbb$read3_code=="K5B00" | ukbb$read3_code=="K5B01" | ukbb$read3_code=="K5B0z")),1,0)

# Identify infertility cases due to anatomical causes based on ICD10 (N971,N972, N973), ICD9 (6283, 6284) and primary care codes (read2: K5B2., K5B20, K5B21, K5B23, K5B2z, K5B3., K5B30, K5B31, K5B3z, K5B4., K5B40, K5B41, K5B5., K5B51; read3: K5B2., K5B20, K5B21, X4099, K5B2z, K5B3., K5B30, K5B31, K5B3z, K5B4., K5B40, K5B41, K5B5., K5B51). 37 cases. 
ukbb$anatomical_infertility<-ifelse((!is.na(ukbb$prim_ICD10) & (ukbb$prim_ICD10=="N971" | ukbb$prim_ICD10=="N972" | ukbb$prim_ICD10=="N973")) | 
(!is.na(ukbb$sec_ICD10) & (ukbb$sec_ICD10=="N971" | ukbb$sec_ICD10=="N972" | ukbb$sec_ICD10=="N973")) |
(!is.na(ukbb$sec_ICD9) & (ukbb$sec_ICD9==6283 | ukbb$sec_ICD9==6284)) |
(!is.na(ukbb$read2) & (ukbb$read2_code=="K5B2." | ukbb$read2_code=="K5B20" | ukbb$read2_code=="K5B21" | ukbb$read2_code=="K5B23" | ukbb$read2_code=="K5B2z" | ukbb$read2_code=="K5B3." | ukbb$read2_code=="K5B30" | ukbb$read2_code=="K5B31" | ukbb$read2_code=="K5B3z" | ukbb$read2_code=="K5B4." | ukbb$read2_code=="K5B40" | ukbb$read2_code=="K5B41" | ukbb$read2_code=="K5B5." | ukbb$read2_code=="K5B51")) | 
(!is.na(ukbb$read3) & (ukbb$read3_code=="K5B2." | ukbb$read3_code=="K5B20" | ukbb$read3_code=="K5B21" | ukbb$read3_code=="X4099" | ukbb$read3_code=="K5B2z" | ukbb$read3_code=="K5B3." | ukbb$read3_code=="K5B30" | ukbb$read3_code=="K5B31" | ukbb$read3_code=="K5B3z" | ukbb$read3_code=="K5B4." | ukbb$read3_code=="K5B40" | ukbb$read3_code=="K5B41" | ukbb$read3_code=="K5B5." | ukbb$read3_code=="K5B51")),1,0)


# Identify idiopathic infertility cases based on code for infertility for other/unknown causes: ICD9: 6288; ICD10: N978; read2: K5By., K5Byz, Kyu9G; read3: K5By., XE0ey, Kyu9G. 
ukbb$idiop_infertility_inclusion<-ifelse((!is.na(ukbb$prim_ICD10) & ukbb$prim_ICD10=="N978") | (!is.na(ukbb$sec_ICD10) & ukbb$sec_ICD10=="N978") | (!is.na(ukbb$sec_ICD9) & ukbb$sec_ICD9==6288) | (!is.na(ukbb$read2) & (ukbb$read2_code=="K5By." | ukbb$read2_code=="K5Byz" | ukbb$read2_code=="Kyu9G"))
| (!is.na(ukbb$read3) & (ukbb$read3_code=="K5By." | ukbb$read3_code=="XE0ey" | ukbb$read3_code=="Kyu9G")),1,0)

# Identify idiopathic infertilitty cases by excluding endometriosis, PCOS, anovulatory infertility and anatomical causes from infertility cases of any cause. 
endo<-read.table("female_ukbb_participants_endometriosis_diagnoses.txt",header=T,stringsAsFactors=F)
endo<-endo[,c(1,16)]
pcos<-read.table("/well/lindgren/UKBIOBANK/laura/pcos/summary_pcos_cases_diagnosed_code_selfrep_symptoms.txt",sep="\t",header=T,stringsAsFactors=F)

ukbb$has_endo<-ifelse(ukbb$eid %in% endo$eid[endo$endometriosis==1],1,0)
ukbb$has_pcos<-ifelse(ukbb$eid %in% pcos$eid,1,0)

ukbb$idiop_infertility_exclusion<-ifelse(ukbb$female_infertility==1 & ukbb$has_endo==0 & ukbb$has_pcos==0 & ukbb$anatomical_infertility==0 & ukbb$anovulatory_infertility==0,1,0) 

# Endometriosis-associated infertility (FinnGen)
ukbb$endometr_infertility<-ifelse(ukbb$has_endo==1 & ukbb$female_infertility==1,1,0)

#Female infertility: cervical, vaginal, other and unspecified. (FinnGen)
ifelse((!is.na(ukbb$prim_ICD10) & (ukbb$prim_ICD10=="N973" | ukbb$prim_ICD10=="N974" | ukbb$prim_ICD10=="N978"| ukbb$prim_ICD10=="N979")) | (!is.na(ukbb$sec_ICD10) & (ukbb$sec_ICD10=="N973" | ukbb$sec_ICD10=="N974" | ukbb$sec_ICD10=="N978"| ukbb$sec_ICD10=="N979")) | (!is.na(ukbb$sec_ICD9) & (ukbb$sec_ICD9==6289 | ukbb$sec_ICD9==6288 | ukbb$sec_ICD9==6284)) | 

# Tubular infertility (FinnGen)



setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(ukbb[,c(1,3,16,17,18,19,22)],"female_ukbb_participants_infertility_cases_5_categories.txt",sep="\t",quote=F,row.names=F)

# Drop non-Europeans from pheno input file
data<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases/female_ukbb_participants_infertility_cases_5_categories.txt",sep="\t",header=T,stringsAsFactors=F)

eur<-read.table("/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt",header=T,sep="\t",stringsAsFactors=F)
eur<-eur[eur$genetically_european==1,]

data<-data[which(data$eid %in% eur$eid),]

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(data,"female_ukbb_participants_infertility_cases_5_categories_europeans_only.txt",sep="\t",quote=F,row.names=F)















