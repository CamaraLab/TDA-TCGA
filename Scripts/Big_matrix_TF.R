suppressWarnings({
  suppressMessages ({
    library(data.table,quietly = T,warn.conflicts = FALSE)
    library(rhdf5,quietly = T,warn.conflicts = FALSE)
    library(getopt,quietly = T,warn.conflicts = FALSE)
  })
  
})

PROJECT_NAME<-"LUADTF"

print ("Reading matrix files:")
TPM.matrix<-fread("LUAD_Full_TPM_matrix.csv",data.table = F,header = T)
rownames(TPM.matrix)<-TPM.matrix[,1]; TPM.matrix<-TPM.matrix[,-1]

mat_non_syn<-fread("tf_mut_total2.csv",data.table = F,header = T)
samples<-mat_non_syn[,1]
mat_non_syn<-mat_non_syn[,-1]
mat_non_syn<-sapply(mat_non_syn,as.numeric)

rownames(mat_non_syn)<-samples
colnames(mat_non_syn)<-paste0(colnames(mat_non_syn),"|0")

clean_samples<-function(matrix) {
  samples<-substring(rownames(matrix),1,15)
  samples_to_keep_1<-which(sapply(samples, function (x) substring(x,14) %in% c("01","10","11","06","07")))
  #Removing duplicated records
  duplicated<-samples[duplicated(samples)]
  samples_to_keep_2<-which(!samples %in% duplicated)
  samples_to_keep<-samples[intersect(samples_to_keep_1,samples_to_keep_2)]
  matrix<-matrix[samples_to_keep,]
  rownames(matrix)<-samples_to_keep
  matrix<-matrix[sort(rownames(matrix)),sort(colnames(matrix))]
  return(as.matrix(matrix))
}





#Removing recurrent tumor and dual samples
print ("Cleaning matrix files:")
TPM.matrix<-clean_samples(TPM.matrix) 

rownames(TPM.matrix)<-substring(rownames(TPM.matrix),1,12)
samples_of_interest<-intersect(rownames(TPM.matrix),rownames(mat_non_syn))

mat_non_syn<-mat_non_syn[samples_of_interest,]
mat_non_syn_bin<-mat_syn<-mat_non_syn
mat_non_syn_bin<-ifelse(mat_non_syn_bin>0,1,0)
mat_syn<-ifelse(mat_syn!=(-1),0,0)

TPM.matrix<-TPM.matrix[samples_of_interest,]

colnames(TPM.matrix)<-paste0("exp_",colnames(TPM.matrix))

print ("Writing Big_matrix file")

BIG.matrix<-cbind(TPM.matrix,mat_non_syn_bin)
write.csv(BIG.matrix,paste0(PROJECT_NAME,"_BIG_matrix.csv"))

print ("Writing h5 files")

h5file<-paste0(PROJECT_NAME,".h5")
if (file.exists(h5file)) {file.remove(h5file)}
h5createFile(h5file)
H5close()
suppressWarnings({
  h5write(mat_non_syn_bin, h5file,"Mutations_Binary",)
  h5write(mat_non_syn,h5file,"Mutations_NS")
  h5write(mat_syn,h5file,"Mutations_S")
  h5write(rownames(mat_non_syn_bin),h5file,"Mutations_Samples")
  h5write(colnames(mat_non_syn_bin),h5file,"Mutations_Genes")
  h5write(TPM.matrix,h5file,"TPM")
  h5write(rownames(TPM.matrix),h5file,"TPM_Samples")
  h5write(colnames(TPM.matrix),h5file,"TPM_Genes")
  
})  

print(h5ls(h5file))
H5close()
