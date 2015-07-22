require(data.table)
require(rhdf5)
setwd("~/TCGA-DATA/GBM")
PROJECT_NAME<-"GBM"

TPM.matrix<-fread("Expression/GBM_Full_TPM_matrix.csv",data.table = F,header = T)
rownames(TPM.matrix)<-TPM.matrix[,1]; TPM.matrix<-TPM.matrix[,-1]

mat_non_syn_bin<-fread("Mutations/GBM_Full_Mutations_binary.csv",data.table = F,header = T)
rownames(mat_non_syn_bin)<-mat_non_syn_bin[,1]; mat_non_syn_bin<-mat_non_syn_bin[,-1]
mat_non_syn<-fread("Mutations/GBM_Full_Mutations_non_synonymous.csv",data.table = F,header = T)
rownames(mat_non_syn)<-mat_non_syn[,1]; mat_non_syn<-mat_non_syn[,-1]
mat_syn<-fread("Mutations/GBM_Full_Mutations_synonymous.csv",data.table = F,header = T)
rownames(mat_syn)<-mat_syn[,1]; mat_syn<-mat_syn[,-1]


TPM.matrix<-clean_samples(TPM.matrix) #Removing recurrent tumor and dual samples
mat_non_syn_bin<-clean_samples(mat_non_syn_bin)
mat_non_syn<-clean_samples(mat_non_syn)
mat_syn<-clean_samples(mat_syn) 


samples_of_interest<-intersect(rownames(TPM.matrix),rownames(mat_non_syn_bin))
TPM.matrix<-TPM.matrix[samples_of_interest,]
mat_non_syn_bin<-mat_non_syn_bin[samples_of_interest,]
mat_non_syn<-mat_non_syn[samples_of_interest,]
mat_syn<-mat_syn[samples_of_interest,]

colnames(mat_non_syn_bin)<-paste0("mut_",colnames(mat_non_syn_bin))
colnames(TPM.matrix)<-paste0("exp_",colnames(TPM.matrix))


BIG.matrix<-cbind(TPM.matrix[samples_of_interest,],mat_non_syn_bin[samples_of_interest,])
write.csv(BIG.matrix,paste0(PROJECT_NAME,"_BIG_matrix.csv"))



h5file<-paste0("Mutations/",PROJECT_NAME,".h5")
h5createFile(h5file)
H5close()
suppressWarnings({
  h5write(mat_non_syn_bin, h5file,"Mutations_Binary")
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
