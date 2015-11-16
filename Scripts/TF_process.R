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


a<-lapply(variant_files,function (x) {
  f<-read.table(x,header=TRUE)[,1:2] 
  f$chrom<-paste0("chr",f$chrom)
  return(f)
})
names(a)<-var_samples


TPM<-read.csv("../../LUAD_CUR/Expression/LUAD_Full_TPM_matrix.csv",row.names = 1,as.is=T)
tpm_samples<-substring(rownames(TPM),1,12)
samples_of_interest<-intersect(tpm_samples,var_samples)

all_chromosomes<-sapply(1:22,function(x) paste0("chr",x))
noncode<-read.table("../NONCODEv5_hg19_linc.bed",as.is=T)[,1:4]
colnames(noncode) <- c('chr','start','end','id')
head(noncode)
noncode<-subset(noncode,noncode$chr %in% all_chromosomes)
dim(noncode)


for (sample in names(a)) {
  print(sample)
  var<-a[[sample]]
  noncode[,sample]<-apply(noncode,1, function (n) {    
    chr<-n[1]
    pos_start<-as.numeric(n[2])
    pos_end<-as.numeric(n[3])
    ans<-in_range(var$position[var$chrom==chr],pos_start,pos_end)
    return(ans)
  })
}

  
noncode_mutation_matrix<-t(noncode[,5:ncol(noncode)])
colnames(noncode_mutation_matrix)<-paste0("mut_",noncode[,4])
noncode_binary_mutation_matrix<-ifelse(noncode_mutation_matrix>0,1,0)
write.csv(noncode_mutation_matrix,"nonocode_mut_total.csv")
write.csv(noncode_binary_mutation_matrix,"nonocode_mut_bin.csv")



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
