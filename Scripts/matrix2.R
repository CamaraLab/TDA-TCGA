# Environemnt settings
#setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
#library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(parallel)
require(getopt,quietly = T)

anno<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv",as.is=T)
rownames(anno)<-anno[,1]


PROJECT_NAME<-"GBM"
wd<-paste0("c:/Users/Udi/Documents/TCGA-DATA/",PROJECT_NAME)
setwd(wd)

#Raw file list.
files<-list.files("./Expression/rsem.genes.results/",full.names=T)
#Creating index file
index.file<-read.delim("./Expression/unc.edu_GBM.IlluminaHiSeq_RNASeqV2.1.4.0.sdrf.txt",as.is=T)
index.cols<-c("Extract.Name","Comment..TCGA.Barcode.")
index<-unique(index.file[,index.cols]) # Remove duplicates
colnames(index)<-c("File","PatientID")
index<-arrange(index,File)

#Creating gene list and file (gene_id)
gene_id_raw<-read.table(files[1],header=T,stringsAsFactors = F)[,"gene_id"] #Reading gene names from first file
gene_id_split<-strsplit(gene_id_raw,"|",fixed=TRUE)
symbol<-sapply(gene_id_split,"[[",1)
id<-sapply(gene_id_split,"[[",2)
old_ids<-intersect(id,anno$Old_ID) #Identifying old ids
id[match(old_ids,id)]<-anno$EntrezID[match(old_ids,anno$Old_ID)] #Replacing old_ids with new ones

symbol<-anno[id,"Symbol"] #Updating symbols avvording to Id's
gene_id<-as.data.frame(cbind(id,symbol))

# Reading scaled expression level into scale.estimates matrix
cl <- makeCluster(detectCores())
clusterExport(cl=cl, varlist=c("files"))
scale.estimates<-parSapply(cl,files,function (x) return(read.table(x,header=T)[,3]))
stopCluster(cl)
scale.estimates<-t(scale.estimates)




#Creating TPM matrix - 
samples<-index$PatientID
TPM.matrix<-as.data.frame(round(log2(1+scale.estimates*10^6),4))
rownames(TPM.matrix)<-samples
colnames(TPM.matrix)<-paste0(gene_id$symbol,"|",gene_id$id)
TPM.matrix<-TPM.matrix[sort(rownames(TPM.matrix)),sort(colnames(TPM.matrix))]
write.csv(TPM.matrix,paste0("Expression/",PROJECT_NAME,"_Full_TPM_matrix.csv"))


