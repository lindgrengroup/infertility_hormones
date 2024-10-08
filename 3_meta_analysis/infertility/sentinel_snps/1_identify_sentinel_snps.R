# setwd("/well/lindgren/laura/projects/infertility/sentinel_snps")
setwd("/well/lindgren/samvida/hormones_infertility/infertility_meta_mvp")

args = commandArgs(trailingOnly=TRUE)

TRAIT=args[1]

if (TRAIT %in% c("female_infertility_analysis1_all", "male_infertility_all")) {
  
  INPUT_FILE<-paste0("filtered/", TRAIT,"_filtered.txt")
  OUTPUT_FILE1<-paste("lead_snps/", TRAIT, "_sentinel_SNPs.txt")
  
  data<-read.table(INPUT_FILE,header=T,stringsAsFactors=F)
  
  # Find sentinel SNPs at genome-wide significance, using distance based clumping 500kb around sentinel SNP
  data_sign<-data[data$PVALUE<=5E-8,]
  
  sentinels<-data.frame(matrix(nrow = 50, ncol = ncol(data)))
  names(sentinels)<-names(data)
  
  for (i in 1:nrow(sentinels)){
    if (nrow(data_sign)>0){
      data_sign<-data_sign[order(data_sign$PVALUE),]
      sentinels[i,]<-data_sign[1,]
      CHR<-data_sign$CHROM[1]
      START_POS<-ifelse(data_sign$GENPOS[1]-500000>0,data_sign$GENPOS[1]-500000,0)
      END_POS<-ifelse(data_sign$GENPOS[1]+500000>0,data_sign$GENPOS[1]+500000,0)
      data_sign$drop<-ifelse(data_sign$CHROM==CHR & data_sign$GENPOS>=START_POS & data_sign$GENPOS<=END_POS,1,0)
      data_sign<-data_sign[data_sign$drop==0,]
      data_sign<-data_sign[,-which(names(data_sign) %in% "drop")]
    }
  }
  
  write.table(sentinels[!is.na(sentinels$ID),],OUTPUT_FILE1,sep="\t",quote=F,row.names=F)
  
  
} else {
  
  INPUT_FILE<-paste(TRAIT,"MA_results_chr_pos.txt",sep="_")
  OUTPUT_FILE1<-paste(TRAIT,"sentinel_SNPs.txt",sep="_")
  # OUTPUT_FILE2<-paste(TRAIT,"sentinel_SNPs_nom_sign.txt")
  
  setwd("/well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process")
  data<-read.table(INPUT_FILE,header=T,stringsAsFactors=F)
  
  
  # Find sentinel SNPs at genome-wide significance, using distance based clumping 500kb around sentinel SNP
  data_sign<-data[data$P.value<=5E-8,]
  
  
  sentinels<-data.frame(matrix(nrow = 50, ncol = ncol(data)))
  names(sentinels)<-names(data)
  
  for (i in 1:nrow(sentinels)){
    if (nrow(data_sign)>0){
      data_sign<-data_sign[order(data_sign$P.value),]
      sentinels[i,]<-data_sign[1,]
      CHR<-data_sign$chr[1]
      START_POS<-ifelse(data_sign$pos[1]-500000>0,data_sign$pos[1]-500000,0)
      END_POS<-ifelse(data_sign$pos[1]+500000>0,data_sign$pos[1]+500000,0)
      data_sign$drop<-ifelse(data_sign$chr==CHR & data_sign$pos>=START_POS & data_sign$pos<=END_POS,1,0)
      data_sign<-data_sign[data_sign$drop==0,]
      data_sign<-data_sign[,-which(names(data_sign) %in% "drop")]
    }
  }
  
  setwd("/well/lindgren/laura/projects/infertility/sentinel_snps")
  write.table(sentinels[!is.na(sentinels$MarkerName),],OUTPUT_FILE1,sep="\t",quote=F,row.names=F)
  
  
  
}


