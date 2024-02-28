# module load R/3.6.2-foss-2019b


setwd("/well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL")

require('coloc')
args = commandArgs(trailingOnly=TRUE)

TRAIT<-args[1]
TISSUE<-args[2]


logfile = paste(TRAIT,paste(TISSUE,"logfile.txt",sep="_"),sep="_")
setwd("/well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL/logs")

sink(logfile)

# Load full gtex gene-variant pairs file
gtexfile<-paste(TISSUE,"allpairs.txt.gz",sep=".")
setwd("/well/lindgren/resources/gtex/v8/GTEx_Analysis_v8_eQTL_all_associations")
gtex<-read.table(gzfile(gtexfile),header=T,stringsAsFactors=F)

# Remove version from gene ID
gtex$gene_id<-gsub("\\..*","",gtex$gene_id)

# Load all genes with TSS within 1Mb of infertility sentinel SNP
setwd("/well/lindgren/laura/projects/infertility/colocalisations")
genes<-read.table("genes_with_TSS_within_1Mb_of_infertility_sentinel_SNP.txt",header=T,sep="\t",stringsAsFactors=F)

# Only keep genes in GTEx file which are also genes data frame. 
gtex<-gtex[which(gtex$gene_id %in% genes$Gene.stable.ID),]

# Add chr, pos, alleles and MarkerName to gtex file
chrpos<-data.frame(do.call('rbind', strsplit(as.character(gtex$variant_id),'_',fixed=TRUE)))
chrpos<-chrpos[,1:4]
names(chrpos)<-c("chr","pos","ref","alt")
chrpos$alt<-as.character(chrpos$alt)
chrpos$ref<-as.character(chrpos$ref)

chrpos$chr<-substring(chrpos$chr,4,nchar(as.character(chrpos$chr)))
chrpos$chr<-as.numeric(chrpos$chr)
chrpos$pos<-as.numeric(as.character(chrpos$pos))
gtex<-cbind(gtex,chrpos)

gtex$alleles<-apply(cbind(gtex$ref, gtex$alt), 1, function(x) paste(sort(x), collapse="_"))
gtex$chrpos<-paste(gtex$chr,gtex$pos,sep=":")
gtex$MarkerName<-paste(gtex$chrpos,gtex$alleles,sep=":")
gtex$N<-gtex$ma_count/(2*gtex$maf)

gtex<-gtex[,c(1,16,13,12,8,9,7,6,17)]
# Load GWAS summary stats for infertility trait
gwasfile<-paste(TRAIT,"MA_results_chr_pos.txt",sep="_")
setwd("/well/lindgren/laura/projects/infertility/meta_analysis/MA_June2023/process")
gwas<-read.table(gwasfile,header=T,stringsAsFactors=F)
gwas<-gwas[,c(1,2,3,4,5,6,7,11,9,10)]

# Only keep overlapping variants in gwas and gtex file, then merge
gwas<-gwas[which(gwas$MarkerName %in% gtex$MarkerName),]
gtex<-gtex[which(gtex$MarkerName %in% gwas$MarkerName),]

names(gwas)<-c("MarkerName","EA_gwas","OA_gwas","EAF_gwas","BETA_gwas" ,"SE_GWAS","P_GWAS","MAF_GWAS"
,"N_CASES","N_CONTROLS")
names(gtex)<-c("gene_id","MarkerName","EA_gtex","OA_gtex" ,"BETA_gtex","SE_gtex", "P_gtex", "MAF_gtex","N_gtex")

data<-merge(gtex,gwas,by="MarkerName",all=T)

# Load sentinel SNPs for infertility trait
sentinelfile<-paste(TRAIT,"sentinel_SNPs.txt",sep="_")
setwd("/well/lindgren/laura/projects/infertility/sentinel_snps")
sentinels<-read.table(sentinelfile,header=T,stringsAsFactors=F)

# Only keep genes with TSS within 1Mb of sentinel variant for this infertility trait
genes$keep2<-0
for (i in 1:nrow(genes)){
	CHR<-genes$Chromosome.scaffold.name[i]
	TSS<-genes$Transcription.start.site..TSS.[i]
	subset<-sentinels[sentinels$chr==CHR & sentinels$pos>TSS-1000000 & sentinels$pos<TSS+1000000,]
	if (nrow(subset)>0){
		genes$keep2[i]<-1
		}
	}

genes<-genes[genes$keep2==1,]
data<-data[which(data$gene_id %in% genes$Gene.stable.ID),]

# Align alleles and effect sizes between gtex and gwas

# Run coloc for all unique ensembl gene IDs in data
setwd("/well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL")
data$EA_gwas<-toupper(data$EA_gwas)
data$OA_gwas<-toupper(data$OA_gwas)

data$aligned<-ifelse(data$EA_gwas==data$EA_gtex & data$OA_gwas==data$OA_gtex,1,
		ifelse(data$EA_gwas==data$OA_gtex & data$OA_gwas==data$EA_gtex,-1,0))

for (i in 1:nrow(data)){
	EA<-data$EA_gwas[i]
	OA<-data$OA_gwas[i]
	BETA<-data$BETA_gwas[i]
	if (data$aligned[i]==-1){
		data$EA_gwas[i]<-OA
		data$OA_gwas[i]<-EA
		data$BETA_gwas[i]<-BETA*(-1)
		}
	}

data<-data[,-which(names(data) %in% "EAF_gwas")]

# Calculate variance
data$varbeta_GWAS<-(data$SE_GWAS)^2
data$varbeta_gtex<-(data$SE_gtex)^2
data<-data[!is.na(data$SE_gtex),]

# Calculate N and case fraction for gWAS
data$N=data$N_CASES+data$N_CONTROLS
data$case_fraction<-data$N_CASES/(data$N_CASES+data$N_CONTROLS)

data<-data[data$MAF_gtex>0.01 & data$MAF_GWAS>0.01,]

# Make data frame to store coloc results
results<-as.data.frame(matrix(ncol=9,nrow=0))
names(results)<-c("nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","gene","Tissue","Trait")


# 1:nrow(filelist)
gene_ids<-unique(data$gene_id)

data<-data[!is.na(data$MAF_gtex) & data$MAF_gtex>0 & data$MAF_gtex<=0.5 & !is.na(data$MAF_GWAS) & data$MAF_GWAS>0 & data$MAF_GWAS<=0.5,]
 
for (i in 1:length(gene_ids)){
	GENE_ID<-gene_ids[i]
	subset<-data[data$gene_id==GENE_ID,]

	# Run coloc
	coloc_df <- coloc.abf(dataset1 = list(snp = subset$MarkerName, pvalues = subset$P_GWAS, beta = subset$BETA_gwas, varbeta = subset$varbeta_GWAS, MAF =subset$MAF_GWAS, s=subset$case_fraction,N = subset$N, type = "cc"),dataset2 = list(snp = subset$MarkerName, pvalues = subset$P_gtex, beta = subset$BETA_gtex, varbeta = subset$varbeta_gtex, MAF =subset$MAF_gtex, N = subset$N_gtex, sdY=1, type = "quant"), p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)

	# Extract and append results  	
  	results[i,1:6]<-coloc_df$summary
	results[i,7]<-GENE_ID
	results[i,8]<-TISSUE
	results[i,9]<-TRAIT
	}

names(results)<-c("nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","gene","Tissue","Trait")

outfile<-paste(TRAIT,paste(TISSUE,"eQTL_GTEx_v8_coloc_results_June2023.txt",sep="_"),sep="_")

setwd("/well/lindgren/laura/projects/infertility/colocalisations/GTEx_eQTL/coloc_results")
write.table(results,outfile,sep="\t",quote=F,row.names=F)

