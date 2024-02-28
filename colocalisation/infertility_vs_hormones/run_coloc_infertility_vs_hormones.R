# module load R/3.6.2-foss-2019b


setwd("/well/lindgren/laura/projects/infertility/colocalisations/sex_hormones")

require('coloc')

args = commandArgs(trailingOnly=TRUE)

INFERTILITY<-args[1]
HORMONE<-args[2]


# Load Infertility GWAS results
infert_infile<-paste(INFERTILITY,"MA_results_chr_pos.txt",sep="_")
setwd("/well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process")
infert<-read.table(infert_infile,header=T,stringsAsFactors=F)
infert<-infert[,c(1,2,3,4,5,6,7,11,9,10)]
infert$Allele1<-toupper(infert$Allele1)
infert$Allele2<-toupper(infert$Allele2)
names(infert)<-c("MarkerName","EA_infert","OA_infert","EAF_infert","BETA_infert",
"SE_infert","P_infert","maf","N_CASES","N_CONTROLS")

# Load Hormone GWAS results
horm_infile<-paste(HORMONE,"1.tbl",sep="_")
setwd("/well/lindgren/laura/projects/infertility/hormones/results")
horm<-read.table(horm_infile,header=T,stringsAsFactors=F)
horm<-horm[,which(names(horm) %in% c("MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","P.value" ,"N"))]
horm$Allele1<-toupper(horm$Allele1)
horm$Allele2<-toupper(horm$Allele2)

horm$alleles<-apply(cbind(horm$Allele1, horm$Allele2), 1, function(x) paste(sort(x), collapse="_"))
horm$MarkerName<-paste(substring(horm$MarkerName,4,nchar(horm$MarkerName)),horm$alleles,sep=":")
horm<-horm[,-which(names(horm) %in% "alleles")]
names(horm)<-c("MarkerName","EA_horm","OA_horm","EAF_horm","BETA_horm","SE_horm","P_horm","N_horm")

# Combine hormone and infertility summary stats in 1 data frame, only keeping variants present in both. 
horm<-horm[which(horm$MarkerName %in% infert$MarkerName),]
infert<-infert[which(infert$MarkerName %in% horm$MarkerName),]

horm<-horm[order(horm$MarkerName),]
infert<-infert[order(infert$MarkerName),]

infert<-merge(infert,horm,by="MarkerName",all=F)
rm(horm)

chrpos<- data.frame(do.call('rbind', strsplit(as.character(infert$MarkerName),':',fixed=TRUE)))
chrpos<-chrpos[,1:2]
names(chrpos)<-c("chr","pos")
chrpos$chr<-as.numeric(as.character(chrpos$chr))
chrpos$pos<-as.numeric(as.character(chrpos$pos))

infert<-cbind(chrpos,infert)

# Load infertility sentinels
setwd("/well/lindgren/laura/projects/infertility/colocalisations/sex_hormones")
sentinels<-read.table("sentinels_infertility_hormones_clumped_100kb.txt",header=T,stringsAsFactors=F)

# Calculate variance
infert$varbeta_infert<-(infert$SE_infert)^2
infert$varbeta_horm<-(infert$SE_horm)^2


# Calculate N and case fraction for infertility
infert$N_infert=infert$N_CASES+infert$N_CONTROLS
infert$case_fraction<-infert$N_CASES/(infert$N_CASES+infert$N_CONTROLS)


# Calculate MAF for homrones
infert$maf_horm<-ifelse(infert$EAF_horm<=0.5,infert$EAF_horm,1-infert$EAF_horm)

# Only keep variants with MAF above 1% in both datasets. 
infert<-infert[infert$maf>0.01 & infert$maf_horm>0.01,]


# Make data frame to store coloc results
results<-as.data.frame(matrix(ncol=10,nrow=0))
names(results)<-c("nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","chr_sentinel","pos_sentinel","infertility","hormone")

for (i in 1:nrow(sentinels)){
	SENTINEL_CHR<-sentinels$chr[i]
	SENTINEL_POS<-sentinels$pos[i]
	subset<-infert[infert$chr==SENTINEL_CHR & infert$pos>SENTINEL_POS-500000 & infert$pos<SENTINEL_POS+500000,]

	# Run coloc
	coloc_df <- coloc.abf(dataset1 = list(snp = subset$MarkerName, pvalues = subset$P_infert, beta = subset$BETA_infert, varbeta = subset$varbeta_infert, MAF =subset$maf, s=subset$case_fraction,N = subset$N, type = "cc"),dataset2 = list(snp = subset$MarkerName, pvalues = subset$P_horm, beta = subset$BETA_horm, varbeta = subset$varbeta_horm, MAF =subset$maf_horm, N = subset$N_horm, sdY=1, type = "quant"), p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)

	# Extract and append results  	
  	results[i,1:6]<-coloc_df$summary
	results[i,7]<-SENTINEL_CHR
	results[i,8]<-SENTINEL_POS
	results[i,9]<-INFERTILITY
	results[i,10]<-HORMONE
	}


# Output results. 
setwd("/well/lindgren/laura/projects/infertility/colocalisations/sex_hormones/results")
outfile<-paste("coloc_results",paste(HORMONE,paste(INFERTILITY,"txt",sep="."),sep="_"),sep="_")
write.table(results,outfile,sep="\t",quote=F,row.names=F)

