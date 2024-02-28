setwd("/well/lindgren/laura/projects/infertility/easyQC")

args = commandArgs(trailingOnly=TRUE)

STUDY<-args[1]
TRAIT<-args[2]
INPUT_FILE<-args[3]
N_CASES<-args[4]
N_CONTROLS<-args[5]

N_CASES<-as.numeric(N_CASES)
N_CONTROLS<-as.numeric(N_CONTROLS)

print(paste("N_CASES",N_CASES,sep=": "))
print(paste("N_CONTROLS",N_CONTROLS,sep=": "))


data<-read.table(INPUT_FILE,header=F,sep="\t",stringsAsFactors=F)
names(data)<-c("chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta",
"af_alt","af_alt_cases","af_alt_controls")

print(paste("Number of SNPs in file pre-QC",nrow(data),sep=":"))


# Drop variants with MAF<1%
# print(paste("Number of SNPs dropped by MAF>1% filter",nrow(data)-nrow(data[data$af_alt<=0.99 & data$af_alt>=0.01,]),sep=":")) 
# data<-data[data$af_alt<=0.99 & data$af_alt>=0.01,]


# Drop variants with extreme or missing beta
print(paste("Number of SNPs dropped by extreme and missing BETA filter",nrow(data)-nrow(data[!is.na(data$beta) & abs(data$beta)<10,]),sep=":"))
data<-data[!is.na(data$beta) & abs(data$beta)<10,]


# Drop variants with extreme, missing or negative SE
print(paste("Number of SNPs by extreme and missing SE filter",nrow(data)-nrow(data[!is.na(data$sebeta) & data$sebeta<10 & data$sebeta>0,]),sep=":"))
data<-data[!is.na(data$sebeta) & data$sebeta<10 & data$sebeta>0,]


# Drop variants with missing LOG10P
print(paste("Number of SNPs by missing LOG10P filter",nrow(data)-nrow(data[!is.na(data$mlogp),]),sep=":"))
data<-data[!is.na(data$mlogp),]

# Drop duplicate chrpos
data$chrpos<-paste(data$chrom,data$pos,sep=":")

print(paste("Number of SNPs dropped by duplicate filter",nrow(data)-nrow(data[duplicated(data$chrpos,fromLast=T)==F & duplicated(data$chrpos,fromLast=F)==F,]),sep=":"))
data<-data[duplicated(data$chrpos,fromLast=T)==F & duplicated(data$chrpos,fromLast=F)==F,]

# Only keep variants with A/G/C/T alleles
data2<-data[which(data$ref %in% c("A","C","G","T")),]
data2<-data2[which(data2$alt %in% c("A","C","G","T")),]
print(paste("Number of SNPs dropped by allele A/C/G/T filter",nrow(data)-nrow(data2),sep=":"))

data<-data2
rm(data2)

# Calculate Neff and MAC
# data$maf<-ifelse(data$af_alt<=0.5,data$af_alt,1-data$af_alt)
# data$Neff<-(4*as.numeric(N_CASES)*as.numeric(N_CONTROLS))/(as.numeric(N_CASES)+as.numeric(N_CONTROLS))
# data$MAC<-2*data$Neff*data$maf
# data$Ncases<-as.numeric(N_CASES)
# data$Ncontrols<-as.numeric(N_CONTROLS)

# Drop variants with MAC<10
# print(paste("Number of SNPs after MAC>=10 filter",nrow(data)-nrow(data[data$MAC>=10,]),sep=":"))
# data<-data[data$MAC>=10,]

# Make markername for MA
data$alleles<-apply(cbind(data$alt, data$ref), 1, function(x) paste(sort(x), collapse="_"))
data$markername<-paste(data$chrpos,data$alleles,sep=":")

# Output QCed data
setwd("/well/lindgren/laura/projects/infertility/easyQC/input_files/FinnGen")
OUTPUT_FILENAME<-paste(STUDY,paste(TRAIT,"noMACfilter_QCed.txt",sep="_"),sep="_")
data<-data[,-which(names(data) %in% c("alleles","chrpos","af_alt_cases","af_alt_controls"))]
write.table(data,OUTPUT_FILENAME,sep="\t",quote=F,row.names=F)


