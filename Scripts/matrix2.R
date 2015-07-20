# Environemnt settings
#setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
#library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(parallel)

anno<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv",as.is=T)
rownames(anno)<-anno[,1]
#anno_old_new<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Anno_old_new.csv",as.is=T)

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
samples<-strtrim(index$PatientID,15)

#Removing non 01,10,11 Samples
samples_to_keep1<-which(sapply(samples, function (x) substring(x,14) %in% c("01","10","11")))
###
#Removing duplicated records
duplicated<-samples[duplicated(samples)]
samples_to_keep2<-which(!samples %in% duplicated)

samples_to_keep<-intersect(samples_to_keep1,samples_to_keep2)


#Creating TPM matrix - 
TPM.matrix<-as.data.frame(round(log2(1+scale.estimates*10^6),4))
TPM.matrix<-TPM.matrix[samples_to_keep,]
rownames(TPM.matrix)<-samples[samples_to_keep]
colnames(TPM.matrix)<-paste0(gene_id$symbol,"|",gene_id$id)
TPM.matrix<-TPM.matrix[sort(rownames(TPM.matrix)),sort(colnames(TPM.matrix))]

write.csv(TPM.matrix,paste0("Expression/",PROJECT_NAME,"_TPM_matrix.csv"))

#Intersecting mut.aut with TPM.matrix
mut<-intersect.mat(TPM.matrix,mut)

#Adding missing samples from TPM matrix to Mutant matrix with zeroes
mut<-rbind(mut,delta.zero.matrix(TPM.matrix,mut))


#Creating Big matrix
TPM.matrix<-TPM.matrix[order(rownames(TPM.matrix)),]
mut<-mut[order(rownames(mut)),]

colnames(TPM.matrix)<-paste0("exp_",colnames(TPM.matrix))
colnames(mut)<-paste0("mut_",colnames(mut))
BIG.matrix<-cbind(TPM.matrix,mut,cnv)

#Splitting cnv matrix to amplification only and deletion only
#In cnv_amp all deletions are zerod, in cnv_del all amplifications are zerod and values are absolute.
cnv_amp<-cnv
cnv_amp[cnv<0]<-0
cnv_del<-cnv
cnv_del[cnv>0]<-0
cnv_del<-abs(cnv_del)
colnames(cnv_amp)<-paste0("amp_",substring(colnames(cnv),4))
colnames(cnv_del)<-paste0("del_",substring(colnames(cnv),4))

#Matrix output to file
#TPM.matrix<-TPM.matrix[mut_curated_samples,]
#cnv_amp<-cnv_amp[mut_curated_samples,]
#cnv_del<-cnv_del[mut_curated_samples,]
#mut<-mut[mut_curated_samples,]
#BIG.matrix<-BIG.matrix[mut_curated_samples,]
write.csv(TPM.matrix,"TPM_matrix_Curated.csv")
write.csv(cnv_amp,"CNV_AMP_matrix_Curated.csv")
write.csv(cnv_del,"CNV_DEL_matrix_Curated.csv")
write.csv(mut,"Mut_matrix_Curated.csv")
write.csv(BIG.matrix,"BIG_matrix_Curated.csv")
