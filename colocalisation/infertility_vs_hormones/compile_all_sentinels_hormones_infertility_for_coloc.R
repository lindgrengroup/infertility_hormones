infert<-read.table("/well/lindgren/laura/projects/infertility/sentinel_snps/sentinel_snps_all_infertility_analyses.txt",sep="\t",header=T,stringsAsFactors=F)
infert<-infert[,-1]

horm<-read.table("/well/lindgren/laura/projects/infertility/hormones/sentinels/all_sentinel_variants_hormones.txt",sep="\t",header=T,stringsAsFactors=F)
horm<-horm[,-1]
names(horm)[1:2]<-c("chr","pos")

sentinels<-merge(infert,horm,all=T)


clumped<-data.frame(matrix(nrow = 1000, ncol = ncol(sentinels)))
names(clumped)<-names(sentinels)

for (i in 1:nrow(clumped)){
	if (nrow(sentinels)>0){
		sentinels<-sentinels[order(sentinels$chr,sentinels$pos),]
		clumped[i,]<-sentinels[1,]
		CHR<-sentinels$chr[1]
		START_POS<-ifelse(sentinels$pos[1]-100000>0,sentinels$pos[1]-50000,0)
		END_POS<-ifelse(sentinels$pos[1]+100000>0,sentinels$pos[1]+50000,0)
		sentinels$drop<-ifelse(sentinels$chr==CHR & sentinels$pos>=START_POS & sentinels$pos<=END_POS,1,0)
		sentinels<-sentinels[sentinels$drop==0,]
		sentinels<-sentinels[,-which(names(sentinels) %in% "drop")]
	}
}

setwd("/well/lindgren/laura/projects/infertility/colocalisations/sex_hormones")
write.table(clumped[!is.na(clumped$chr),],"sentinels_infertility_hormones_clumped_100kb.txt",sep="\t",quote=F,row.names=F)


