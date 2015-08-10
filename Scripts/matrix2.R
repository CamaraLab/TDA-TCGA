# Environemnt settings
#setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
#library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(parallel)
require(getopt,quietly = T)


PROJECT_NAME<-"LUSC"

anno<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv",as.is=T)
rownames(anno)<-anno[,1]
anno_old_new<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/anno_old_new.csv",as.is=T,row.names=1)


wd<-paste0("c:/Users/Udi/Documents/TCGA-DATA/",PROJECT_NAME)
setwd(wd)

#Raw file list.
files<-list.files("./Expression/rsem",full.names=T)
#Creating index file
index.file<-read.delim("./Expression/unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.1.11.0.sdrf.txt",as.is=T)
index.cols<-c("Extract.Name","Comment..TCGA.Barcode.")
index<-unique(index.file[,index.cols]) # Remove duplicates
colnames(index)<-c("File","PatientID")
index<-arrange(index,File)

#Creating gene list and file (gene_id)
gene_id_raw<-read.table(files[1],header=T,stringsAsFactors = F)[,"gene_id"] #Reading gene names from first file
gene_id_split<-strsplit(gene_id_raw,"|",fixed=TRUE)
symbol<-sapply(gene_id_split,"[[",1)
id<-as.numeric(sapply(gene_id_split,"[[",2))


#Replacing old ids witht new
ids_to_replace_indices<-id %in% anno_old_new$Old_ID
ids_to_replace<-id[ids_to_replace_indices]
new_id_indices<-match(ids_to_replace,anno_old_new$Old_ID)
id[ids_to_replace_indices]<-anno_old_new$New_ID[new_id_indices]

symbol<-anno[as.character(id),"Symbol"] #Updating symbols avvording to Id's
gene_id<-as.data.frame(cbind(id,symbol),stringsAsFactors=F)

#Bring back withdrawn entrez ids
withdrawn_genes_indices<-as.numeric(rownames(gene_id[is.na(gene_id$id),]))
original_id<-sapply(gene_id_split[withdrawn_genes_indices],"[[",2)
gene_id[withdrawn_genes_indices,"id"]<-original_id
gene_id[withdrawn_genes_indices,"symbol"]<-rep("withdrawn",length(original_id))


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
TPM.matrix<-round(TPM.matrix,5)
write.csv(TPM.matrix,paste0("Expression/",PROJECT_NAME,"_Full_TPM_matrix.csv"))


