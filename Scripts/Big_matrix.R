
udiread<-function(file1) {
  x<-fread(file1,data.table=FALSE)
  row.names(x)<-x[,1]
  x<-x[,-1]
}

suppressWarnings({
  suppressMessages ({
    library(data.table,quietly = T,warn.conflicts = FALSE)
    library(rhdf5,quietly = T,warn.conflicts = FALSE)
    library(getopt,quietly = T,warn.conflicts = FALSE)
    require("data.table")
  })
  
})

spec = matrix(c(
  "expression", "e",1, "character",
  "binary", "b",1, "character",
  "syn", "s",1, "character",
  "non_syn", "n",1, "character",
  "cnv","c",2,"character",
  "project", "p", 1, "character",
  "BIG","B",2,"character"
  
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode
print(names(arg))
print (arg)

#setwd("../../../../Udi/SkyDrive/TCGA_CURATED/LAML/")
#arg$expression<-"Expression/LAML_Full_TPM_matrix.csv"
#arg$binary<-"Mutations/LAML_Full_Mutations_binary.csv"
#arg$syn<-"Mutations/LAML_Full_Mutations_synonymous.csv"
#arg$non_syn<-"Mutations/LAML_Full_Mutations_non_synonymous.csv"
#arg$project<-"LAML"

PROJECT_NAME<-arg$project
#wd<-paste0("~/TCGA-DATA/",PROJECT_NAME)
#setwd(wd)
#TPM.matrix<-fread("Expression/GBM_Full_TPM_matrix.csv",data.table = F,header = T)
#TPM.matrix<-fread(paste0("Expression/",PROJECT_NAME,"_Full_TPM_matrix.csv"),data.table = F,header = T)

print ("Reading matrix files:")
if (!("BIG" %in% names(arg))) { #If big matrix is not supplied
  
  
  TPM.matrix<-fread(arg$expression,data.table = F,header = T)
  rownames(TPM.matrix)<-TPM.matrix[,1]; TPM.matrix<-TPM.matrix[,-1]
  
  #mat_non_syn_bin<-fread(paste0("Mutations/",PROJECT_NAME,"_Full_Mutations_binary.csv"),data.table = F,header = T)
  mat_non_syn_bin<-fread(arg$binary,data.table = F,header = T)
  rownames(mat_non_syn_bin)<-mat_non_syn_bin[,1]; mat_non_syn_bin<-mat_non_syn_bin[,-1]
  #mat_non_syn<-fread(paste0("Mutations/",PROJECT_NAME,"_Full_Mutations_non_synonymous.csv"),data.table = F,header = T)
  mat_non_syn<-fread(arg$non_syn,data.table = F,header = T)
  rownames(mat_non_syn)<-mat_non_syn[,1]; mat_non_syn<-mat_non_syn[,-1]
  #mat_syn<-fread(paste0("Mutations/",PROJECT_NAME,"_Full_Mutations_synonymous.csv"),data.table = F,header = T)
  mat_syn<-fread(arg$syn,data.table = F,header = T)
  rownames(mat_syn)<-mat_syn[,1]; mat_syn<-mat_syn[,-1]
  
  
  
  clean_samples<-function(matrix) {
    samples<-substring(rownames(matrix),1,15)
    samples_to_keep_1<-which(sapply(samples, function (x) substring(x,14) %in% c("01","10","11","06","07","03")))
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
  mat_non_syn_bin<-clean_samples(mat_non_syn_bin)
  mat_non_syn<-clean_samples(mat_non_syn)
  mat_syn<-clean_samples(mat_syn) 
  
  #Inersecting samples
  print ("Intersecting samples")
  samples_of_interest<-intersect(rownames(TPM.matrix),rownames(mat_non_syn_bin))
  TPM.matrix<-TPM.matrix[samples_of_interest,]
  mat_non_syn_bin<-mat_non_syn_bin[samples_of_interest,]
  mat_non_syn<-mat_non_syn[samples_of_interest,]
  mat_syn<-mat_syn[samples_of_interest,]
  
  colnames(mat_non_syn_bin)<-paste0("mut_",colnames(mat_non_syn_bin))
  colnames(TPM.matrix)<-paste0("exp_",colnames(TPM.matrix))
  
  
  BIG.matrix<-cbind(TPM.matrix[samples_of_interest,],mat_non_syn_bin[samples_of_interest,])
  filename<-paste0(PROJECT_NAME,"_BIG_matrix.csv")
  
} else { # BIG matrix supplied
  BIG.matrix<-udiread(arg$BIG)
}


 
if ("cnv" %in% names(arg)) {

  print ("Incorporating CNV information")  
  y<-udiread(arg$cnv)
  y<-t(y)
  
  del<-y
  del[del>0]<-0
  del<-del*(-1)
  colnames(del)<-paste0("del_",colnames(del))
  del<-as.data.frame(del)
  
  amp<-y
  amp[amp<0]<-0
  colnames(amp)<-paste0("amp_",colnames(amp))
  amp<-as.data.frame(amp)
  BIG.matrix<-cbind(BIG.matrix,amp[rownames(BIG.matrix),],del[rownames(BIG.matrix),])
  filename<-paste0(PROJECT_NAME,"_BIG_matrix_CNV.csv")
}



print ("Writing Big_matrix file")

#dim(BIG.matrix)
#write.csv(BIG.matrix,paste0(PROJECT_NAME,"_BIG_matrix.csv"))

write.csv(as.matrix(BIG.matrix),filename)



if (!("BIG" %in% names(arg))) {
  print ("Writing h5 files")
  
  h5file<-paste0(PROJECT_NAME,".h5")
  if (file.exists(h5file)) {file.remove(h5file)}
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
  
}
