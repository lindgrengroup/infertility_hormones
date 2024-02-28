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


data<-read.table(INPUT_FILE,header=T,sep=",",stringsAsFactors=F)


print(paste("Number of SNPs in file pre-QC",nrow(data),sep=":"))


# Drop variants with MAF<1%
# print(paste("Number of SNPs dropped by MAF>1% filter",nrow(data)-nrow(data[data$af_alt<=0.99 & data$af_alt>=0.01,]),sep=":")) 
# data<-data[data$af_alt<=0.99 & data$af_alt>=0.01,]


# Drop variants with extreme or missing beta
print(paste("Number of SNPs dropped by extreme and missing BETA filter",nrow(data)-nrow(data[!is.na(data$BETA) & abs(data$BETA)<10,]),sep=":"))
data<-data[!is.na(data$BETA) & abs(data$BETA)<10,]


# Drop variants with extreme, missing or negative SE
print(paste("Number of SNPs by extreme and missing SE filter",nrow(data)-nrow(data[!is.na(data$SE) & data$SE<10 & data$SE>0,]),sep=":"))
data<-data[!is.na(data$SE) & data$SE<10 & data$SE>0,]


# Drop variants with missing LOG10P
print(paste("Number of SNPs by missing p-value filter",nrow(data)-nrow(data[!is.na(data$p.value),]),sep=":"))
data<-data[!is.na(data$p.value),]

# Drop duplicate chrpos
data$chrpos<-paste(data$CHR,data$POS,sep=":")

print(paste("Number of SNPs dropped by duplicate filter",nrow(data)-nrow(data[duplicated(data$chrpos,fromLast=T)==F & duplicated(data$chrpos,fromLast=F)==F,]),sep=":"))
data<-data[duplicated(data$chrpos,fromLast=T)==F & duplicated(data$chrpos,fromLast=F)==F,]

# Only keep variants with A/G/C/T alleles
data2<-data[which(data$Allele1 %in% c("A","C","G","T")),]
data2<-data2[which(data2$Allele2 %in% c("A","C","G","T")),]
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
data$alleles<-apply(cbind(data$Allele1, data$Allele2), 1, function(x) paste(sort(x), collapse="_"))
data$markername<-paste(data$chrpos,data$alleles,sep=":")

# Output QCed data
setwd("/well/lindgren/laura/projects/infertility/easyQC/input_files/Danes")
OUTPUT_FILENAME<-paste(STUDY,paste(TRAIT,"QCed.txt",sep="_"),sep="_")
data<-data[,-which(names(data) %in% c("alleles","chrpos","AF_case","AF_ctrl","N_case","N_ctrl","N_case_hom","N_case_het","N_ctrl_hom","N_ctrl_het"))]
write.table(data,OUTPUT_FILENAME,sep="\t",quote=F,row.names=F)


