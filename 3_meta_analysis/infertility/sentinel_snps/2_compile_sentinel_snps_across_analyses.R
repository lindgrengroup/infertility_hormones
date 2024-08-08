
setwd("/well/lindgren/laura/projects/infertility/sentinel_snps")
fi1_all<-read.table("female_infertility_analysis1_all_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
fi1_eur<-read.table("female_infertility_analysis1_eur_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
fi2_all<-read.table("female_infertility_analysis2_all_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
# fi2_eur<-read.table("female_infertility_analysis2_eur_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
fi3_all<-read.table("female_infertility_analysis3_all_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
fi3_eur<-read.table("female_infertility_analysis3_eur_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
fi4_all<-read.table("female_infertility_analysis4_all_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
fi4_eur<-read.table("female_infertility_analysis4_eur_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
fi5_eur<-read.table("female_infertility_analysis5_eur_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
mi_all<-read.table("male_infertility_all_sentinel_SNPs.txt",header=T,stringsAsFactors=F)
mi_eur<-read.table("male_infertility_eur_sentinel_SNPs.txt",header=T,stringsAsFactors=F)

fi1_all$fi1_all<-1
fi1_eur$fi1_eur<-1
fi2_all$fi2_all<-1
# fi2_eur$fi2_eur<-1
fi3_all$fi3_all<-1
fi3_eur$fi3_eur<-1
fi4_all$fi4_all<-1
fi4_eur$fi4_eur<-1
# fi5_all$fi5_all<-1
fi5_eur$fi5_eur<-1
mi_all$mi_all<-1
mi_eur$mi_eur<-1

fi1_all<-fi1_all[,c(1,12,13,14)]
fi1_eur<-fi1_eur[,c(1,12,13,14)]
fi2_all<-fi2_all[,c(1,12,13,14)]
#fi2_eur<-fi2_eur[,c(1,12,13,14)]
fi3_all<-fi3_all[,c(1,12,13,14)]
fi3_eur<-fi3_eur[,c(1,12,13,14)]
fi4_all<-fi4_all[,c(1,12,13,14)]
fi4_eur<-fi4_eur[,c(1,12,13,14)]
#fi5_all<-fi5_all[,c(1,12,13,14)]
fi5_eur<-fi5_eur[,c(1,12,13,14)]
mi_all<-mi_all[,c(1,12,13,14)]
mi_eur<-mi_eur[,c(1,12,13,14)]

sentinels<-merge(fi1_all,fi1_eur,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,fi2_all,by=c("MarkerName","chr","pos"),all=T)
# sentinels<-merge(sentinels,fi2_eur,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,fi3_eur,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,fi3_all,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,fi4_eur,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,fi4_all,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,fi5_eur,by=c("MarkerName","chr","pos"),all=T)
# sentinels<-merge(sentinels,fi5_all,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,mi_all,by=c("MarkerName","chr","pos"),all=T)
sentinels<-merge(sentinels,mi_eur,by=c("MarkerName","chr","pos"),all=T)

write.table(sentinels,"sentinel_snps_all_infertility_analyses.txt",sep="\t",quote=F,row.names=F)

