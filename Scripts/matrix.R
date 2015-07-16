# Environemnt settings
#setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
library(org.Hs.eg.db)
library(plyr)
library(stringr)
library(parallel)

anno<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv",as.is=T)
anno_old_id<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/entrez_oldid_conversion.csv",as.is=T)

#Raw file list.
files<-list.files("./rsem.genes.results/",full.names=T)
#Creating index file
index.file<-read.delim("unc.edu_GBM.IlluminaHiSeq_RNASeqV2.1.4.0.sdrf.txt")
index.cols<-c("Extract.Name","Comment..TCGA.Barcode.")
index<-unique(index.file[,index.cols]) # Remove duplicates
colnames(index)<-c("File","PatientID")
index<-arrange(index,File)

#Creating gcgene list and file (gene_id)
gene_id_raw<-read.table(files[1],header=T,stringsAsFactors=F)[,"gene_id"] #Reading gene names from first file
gene_id<-as.data.frame(fix_symbols(gene_id_raw),stringsAsFactors = F)
write.csv(gene_id,"gene_table.csv",row.names=FALSE)
gene_id<-gene_id_cur<-read.csv("gene_id_cur.csv")

#Replacing old_entrez_id
unknown_id_index<-which(gene_id$symbol=="?")
unknown_id<-gene_id$id[unknown_id_index]
new_id<-anno_old_id$New_ID[match(unknown_id,anno_old_id$Old_ID)]
new_symbol<-anno$Symbol[match(new_id,anno$EntrezID)]
gene_id$id[unknown_id_index]<-new_id
gene_id$symbol[unknown_id_index]<-new_symbol
# Reading scaled expression level into scale.estimates matrix
cl <- makeCluster(detectCores())
clusterExport(cl=cl, varlist=c("LUAD.files"))
scale.estimates<-parSapply(cl,LUAD.files,function (x) return(read.table(x,header=T)[,3]))
stopCluster(cl)
scale.estimates<-t(scale.estimates)
rownames(scale.estimates)<-strtrim(index$PatientID,15)
colnames(scale.estimates)<-gene_id_cur[,"Symbol"]

#Creating TPM matrix - 
TPM.matrix<-as.data.frame(round(log2(1+scale.estimates*10^6),4))
TPM.matrix<-TPM.matrix[order(rownames(TPM.matrix)),]
#colnames(TPM.matrix)<-paste0("exp_",gene_id[,"Symbol"])
colnames(TPM.matrix)<-gene_id[,"Symbol"]

#Read mutation data
#mut<-t(read.table("./MUTATIONS/UCSC/genomicMatrix.automated_fixed.txt",stringsAsFactors=FALSE,comment.char="",header=T,row.names=1))
mut<-t(read.table("./MUTATIONS/UCSC/genomicMatrix.curated",stringsAsFactors=FALSE,comment.char="",header=T))
colnames(mut)<-mut[1,]
mut<-mut[-1,]
mut<-fix_patient_id(mut)
mut_curated_samples<-rownames(mut)
#Intersecting mut.aut with TPM.matrix
mut<-intersect.mat(TPM.matrix,mut)

#Adding missing samples from TPM matrix to Mutant matrix with zeroes
mut<-rbind(mut,delta.zero.matrix(TPM.matrix,mut))


#Read CNV Data
cnv<-t(read.table("./CNV/genomicMatrix.Gistic2",stringsAsFactors=FALSE,header=T,row.names=1,sep="\t"))
cnv<-fix_patient_id(cnv)

#Intersecting CNV with TPM.matrix
cnv<-intersect.mat(TPM.matrix,cnv)
#Adding missing samples from TPM matrix to Mutant matrix with zeroes
cnv<-rbind(cnv,delta.zero.matrix(TPM.matrix,cnv))


# Read clinical file
clinical<-as.matrix(read.delim("./clinical/clinical/Biotab/nationwidechildrens.org_clinical_patient_luad.txt",stringsAsFactors=FALSE,row.names=2))
clinical<-as.matrix(clinical[3:nrow(clinical),"death_days_to"])
#clinical<-arrange(clinical,bcr_patient_barcode)
#rownames(clinical)<-clinical[,1]
#clinical<-clinical[,-1]

dup_TPM<-duplicated(substr(rownames(TPM.matrix),1,12))
dup_Clinical<-duplicated(rownames(clinical))
not_dup<-rownames(TPM.matrix)[!dup]
a<-intersect(substr(not_dup,1,12),rownames(clinical))
clinical<-as.matrix(clinical[a,])
#Intersecting clinical with TPM.matrix

i<-pmatch(rownames(clinical),rownames(TPM.matrix))
j<-i[!is.na(i)]

clinical<-intersect.mat(TPM.matrix,as.matrix(clinical))
#Adding missing samples from TPM matrix to clinical matrix with zeroes
clinical<-rbind(clinical,delta.zero.matrix(TPM.matrix,clinical))



#Creating Big matrix
TPM.matrix<-TPM.matrix[order(rownames(TPM.matrix)),]
cnv<-cnv[order(rownames(cnv)),]
mut<-mut[order(rownames(mut)),]

colnames(TPM.matrix)<-paste0("exp_",colnames(TPM.matrix))
colnames(mut)<-paste0("mut_",colnames(mut))
colnames(cnv)<-paste0("cnv_",colnames(cnv))
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
