
setwd("/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files")

poi<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases/POI_eur_regenie_europeans_only.txt",header=T,stringsAsFactors=F)

poi_ext<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases/POI_extended_eur_regenie_europeans_only.txt",header=T,stringsAsFactors=F)
poi_ext<-poi_ext[,c(1,2,3)]

ihh<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases/ihh_eur_regenie_europeans_only.txt",header=T,stringsAsFactors=F)

fi<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases/female_ukbb_participants_infertility_cases_5_categories_europeans_only.txt",header=T,stringsAsFactors=F)
fi<-fi[,-which(names(fi) %in% c("age_at_recruitment"))]

mi<-read.table("/well/lindgren/UKBIOBANK/laura/infertility/compile_cases/male_infertility_traits_eur_regenie_europeans_only.txt",header=T,stringsAsFactors=F)

covars<-read.table("covar_file_regenie.txt",header=T,stringsAsFactors=F)

data<-merge(covars,ihh,by=c("FID","IID"),all.x=F,all.y=F)
data$IHH_female<-ifelse(data$sex==2,data$IHH,NA)
data$IHH_male<-ifelse(data$sex==1,data$IHH,NA)

data<-merge(data,poi,by=c("FID","IID"),all.x=T,all.y=F)
data<-merge(data,poi_ext,by=c("FID","IID"),all.x=T,all.y=F)
data<-merge(data,fi,by.x="IID",by.y="eid",all.x=T,all.y=F)
data<-merge(data,mi,by=c("FID","IID"),all.x=T,all.y=F)

covar<-data[,seq(1,29,1)]
pheno<-data[,c(1,2,seq(30,40,1))]

setwd("/well/lindgren/UKBIOBANK/laura/infertility/REGENIE/input_files")
write.table(covar,"regenie_covar_input_file_all_traits_EUR_ONLY.txt",sep="\t",quote=F,row.names=F)
write.table(pheno,"regenie_pheno_input_file_all_traits_EUR_ONLY.txt",sep="\t",quote=F,row.names=F)

