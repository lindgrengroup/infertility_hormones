setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")



cases<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/POI_eur_regenie.txt",sep="\t",header=T,stringsAsFactors=F)
aam<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/self_reported/self_reported_age_at_menopause_3581_0_1_2.txt",sep=",",header=T,stringsAsFactors=F)

aam$f.3581.0.0<-ifelse(!is.na(aam$f.3581.0.0) & aam$f.3581.0.0<0,NA,aam$f.3581.0.0)
aam$f.3581.1.0<-ifelse(!is.na(aam$f.3581.1.0) & aam$f.3581.1.0<0,NA,aam$f.3581.1.0)
aam$f.3581.2.0<-ifelse(!is.na(aam$f.3581.2.0) & aam$f.3581.2.0<0,NA,aam$f.3581.2.0)

aam<-aam[!is.na(aam$f.3581.0.0) | !is.na(aam$f.3581.1.0) | !is.na(aam$f.3581.2.0),]


aam$age_at_menopause<-apply(aam[, 2:4], 1, max,na.rm=T)

data<-merge(cases,aam[,c(1,5)],by.x="IID",by.y="f.eid",all.x=T,all.y=F)
names(data)[which(names(data) %in% "POI")]<-"POI_by_code"

data$POI_by_aam<-ifelse(!is.na(data$age_at_menopause) & data$age_at_menopause>16 & data$age_at_menopause<40 ,1,0)

data$POI_code_and_aam<-ifelse(data$POI_by_aam==1 & data$POI_by_code==1,1,0)

data$POI_code_or_aam<-ifelse(data$POI_by_aam==1 | data$POI_by_code==1,1,0)
data<-data[,-which(names(data) %in% "age_at_menopause")]

data<-data[,c(2,1,4,5,6)]
write.table(data,"/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/POI_extended_eur_regenie.txt",sep="\t",quote=F,row.names=F)


# Drop non-Europeans from pheno input file
data<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files/POI_extended_eur_regenie.txt",sep="\t",header=T,stringsAsFactors=F)

eur<-read.table("/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt",header=T,sep="\t",stringsAsFactors=F)
eur<-eur[eur$genetically_european==1,]

data<-data[which(data$IID %in% eur$eid),]

setwd("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases")
write.table(data,"POI_extended_eur_regenie_europeans_only.txt",sep="\t",quote=F,row.names=F)




