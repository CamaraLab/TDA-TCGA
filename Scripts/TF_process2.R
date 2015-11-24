#setwd("LUAD_WG_HC_filtered/")
library(getopt,quietly = T,warn.conflicts = FALSE)
require("data.table")
require("parallel")
spec = matrix(c(
  "bed_file", "b",1,"character",
  "chunk", "k",2,"numeric",
  "cores", "q",2,"numeric"
), byrow=TRUE, ncol=4)



arg<-getopt(spec) #Conmment this line for debug mode

if ( is.null(arg$chunk ) ) {arg$chunk= 10000}
if ( is.null(arg$cores ) ) {arg$cores= 4}


in_range<-function(value,start,end) {
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


print ("Reading bed file")
all_chromosomes<-sapply(1:22,function(x) paste0("chr",x))
tf<-read.table(arg$bed_file,as.is=T)[,1:4]
#tf<-fread(arg$bed_file,as.is=T)[,1:4]
#tf<-read.table("NONCODEv5_hg19_linc.bed",as.is=T)[,1:4]
#colnames(noncode) <- c('chr','start','end','id')
colnames(tf) <- c('chr','start','end','id')
#tf<-as.matrix(tf)
tf[,"chr"]<-as.character(tf[,"chr"])
tf[,"start"]<-as.numeric(tf[,"start"])
tf[,"end"]<-as.numeric(tf[,"end"])


q<-seq_along(1:nrow(tf)) #Subsetting columns
split.q<-split(q,ceiling(seq_along(q)/arg$chunk))


cl <- makeCluster(arg$cores)
#cl <- makeCluster(2)

print ("Creating mutation matrix")
for (sample in names(a))  {
  print(Sys.time())
  print(sample)
  var<-a[[sample]]
  varlist=c("tf","in_range","var")
  clusterExport(cl=cl, varlist=varlist,envir=environment())
  
  y<-parLapply(cl,split.q,function (x) {
    tf1<-tf[x,]
    mut_column<-apply(tf1,1, function (n) {    
      chr<-n[1]
      pos_start<-as.numeric(n[2])
      pos_end<-as.numeric(n[3])
      ans<-in_range(var$position[var$chrom==chr],pos_start,pos_end)
      return(ans)
    })
    #tf1<-cbind(tf1,mut_column)
    #colnames(tf1)[ncol(tf1)]<-sample
    return(mut_column)
  })
  #mut_column<-apply(tf,1, function (n) {    
  #  chr<-n[1]
  #  pos_start<-as.numeric(n[2])
  #  pos_end<-as.numeric(n[3])
  #  ans<-in_range(var$position[var$chrom==chr],pos_start,pos_end)
  #  return(ans)
  #})
    a1<-NULL
  for (i in 1:length(y)) {
    a1<-c(a1,y[[i]])
  }
  tf<-cbind(tf,a1)
  colnames(tf)[ncol(tf)]<-sample
  gc()
  

  print(Sys.time())
}





tf_mutation_matrix<-t(tf[,5:ncol(tf)])
colnames(tf_mutation_matrix)<-paste0("mut_",tf[,4])
tf_binary_mutation_matrix<-ifelse(tf_mutation_matrix>0,1,0)
write.csv(tf_mutation_matrix,"tf_mut_total2.csv")
write.csv(tf_binary_mutation_matrix,"tf_mut_bin2.csv")



#tf_mutation_matrix<-tf[,5:ncol(tf)]
#rownames(tf_mutation_matrix)<-tf[,4]
#tf_mutation_matrix<-as.data.frame(tf_mutation_matrix,stringsAsFactors = FALSE)
#for (i in 1:ncol(tf_mutation_matrix)) {
#  tf_mutation_matrix[,i]<-as.numeric(tf_mutation_matrix[,i])
#}
##tf_mutation_matrix<-t(tf_mutation_matrix)
#tf_binary_mutation_matrix<-ifelse(tf_mutation_matrix>0,1,0)

#tf_mutation_matrix<-as.matrix(tf_mutation_matrix)
#tf_binary_mutation_matrix<-as.matrix(tf_binary_mutation_matrix)
##write.csv(tf_mutation_matrix,"tf_mut_total2.csv")
#write.csv(tf_binary_mutation_matrix,"tf_mut_bin2.csv")

