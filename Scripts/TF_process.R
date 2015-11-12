#require("GenomicRanges")
#require("IRanges")
#library(GenomicFeatures)

in_range<-function(value,start,end) {
  #x<-(value %in% seq(start,end)) 
  #return(sum(x)) #Returns number of mutations in range
  sum( start <= value & value <= end )
}

setwd("C:/Users/Udi/SkyDrive/TCGA_CURATED/WHOLE_GENOME/LUAD_WG_HC_filtered")
variant_files<-list.files(pattern=".snp.Somatic")
var_samples<-substring(variant_files,1,12)
var_samples<-gsub(pattern = "_",replacement = "-",x=var_samples)

TPM<-read.csv("../../LUAD_CUR/Expression/LUAD_Full_TPM_matrix.csv",row.names = 1,as.is=T)
tpm_samples<-substring(rownames(TPM),1,12)
samples<-intersect(tpm_samples,var_samples)



a<-read.table("TCGA_05_4249.snp.Somatic.hc.filtered",header=TRUE)
noncode<-read.table("../NONCODEv5_hg19_linc.bed")
noncode<-noncode[,1:4]
colnames(noncode) <- c('chr','start','end','id')

x<-NULL
for (i in 1:nrow(noncode)) {
  x[i]<-in_range(a$position,noncode[i,"start"],noncode[i,"end"])
}

#tf<-read.table("../TF_ENCODE_narrowPeak_hg19.bed")
#colnames(tf) <- c('chr','start','end','id')
tf1<-as.data.table(tf)
head(tf)

x<-NULL

for (i in 1:nrow(tf)) {
  x[i]<-in_range(a$position,tf[i,"start"],tf[i,"end"])
}


bed <- with(tf, GRanges(chr, IRanges(start, end),id=id))

View(patients)
