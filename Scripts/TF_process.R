library(getopt,quietly = T,warn.conflicts = FALSE)

spec = matrix(c(
  "bed_file", "b",1,"character"
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode





in_range<-function(value,start,end) {
  #x<-(value %in% seq(start,end)) 
  #return(sum(x)) #Returns number of mutations in range
  sum( start <= value & value <= end )
}

#setwd("C:/Users/Udi/SkyDrive/TCGA_CURATED/WHOLE_GENOME/LUAD_WG_HC_filtered")
variant_files<-list.files(pattern=".snp.Somatic")
var_samples<-substring(variant_files,1,12)
var_samples<-gsub(pattern = "_",replacement = "-",x=var_samples)


print ("Reading variant files")
a<-lapply(variant_files,function (x) {
  f<-read.table(x,header=TRUE)[,1:2] 
  f$chrom<-paste0("chr",f$chrom)
  return(f)
})
names(a)<-var_samples


#TPM<-read.csv("../../LUAD_CUR/Expression/LUAD_Full_TPM_matrix.csv",row.names = 1,as.is=T)
#tpm_samples<-substring(rownames(TPM),1,12)
#samples_of_interest<-intersect(tpm_samples,var_samples)


print ("Reading bed file")
all_chromosomes<-sapply(1:22,function(x) paste0("chr",x))
#tf<-read.table("../TF_ENCODE_narrowPeak_hg19.bed",as.is=T)[,1:4]
tf<-read.table(arg$bed_file,as.is=T)[,1:4]
#noncode<-read.table("../NONCODEv5_hg19_linc.bed",as.is=T)[,1:4]
#colnames(noncode) <- c('chr','start','end','id')
colnames(tf) <- c('chr','start','end','id')

print ("Creating mutation matrix")
for (sample in names(a))  {
  print(Sys.time())
  print(sample)
  var<-a[[sample]]
  tf[,sample]<-apply(tf,1, function (n) {    
    chr<-n[1]
    pos_start<-as.numeric(n[2])
    pos_end<-as.numeric(n[3])
    ans<-in_range(var$position[var$chrom==chr],pos_start,pos_end)
    return(ans)
  })
  print(Sys.time())
}

  
tf_mutation_matrix<-t(tf[,5:ncol(tf)])
colnames(tf_mutation_matrix)<-paste0("mut_",tf[,4])
tf_binary_mutation_matrix<-ifelse(tf_mutation_matrix>0,1,0)
write.csv(tf_mutation_matrix,"tf_mut_total.csv")
write.csv(tf_binary_mutation_matrix,"tf_mut_bin.csv")

