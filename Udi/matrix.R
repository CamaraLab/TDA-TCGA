# Environemnt settings
setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
library(org.Hs.eg.db)
library(plyr)
library(stringr) 

#Raw file list.
LUAD.files<-list.files("./rsem.genes.results/",full.names=T)
#Creating index file
index.file<-read.delim("unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.1.15.0.sdrf.txt")
index.cols<-c("Extract.Name","Comment..TCGA.Barcode.")
index<-unique(index.file[,index.cols]) # Remove duplicates
colnames(index)<-c("File","PatientID")
index<-arrange(index,File)

#Creating gcgene list and file (gene_id)
gene_id_raw<-read.table(LUAD.files[1],header=T,stringsAsFactors=F)[,"gene_id"]
gene_id<-fix_symbols(gene_id_raw)
write.csv(gene_id,"gene_table.csv",row.names=FALSE)
gene_id<-gene_id_cur<-read.csv("gene_id_cur.csv")


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
mut<-t(read.table("./MUTATIONS/UCSC/genomicMatrix.automated_fixed.txt",stringsAsFactors=FALSE,comment.char="",header=T,row.names=1))
mut<-fix_patient_id(mut)
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

#Matrix output to file
write.csv(TPM.matrix,"TPM.matrix.csv")
write.csv(cnv,"CNV.matrix.csv")
write.csv(mut,"Mut.matrix.csv")
write.csv(BIG.matrix,"BIG.matrix3.csv")



########################
# Many bulshit down from here#
#####################33









#cnv[]<-lapply(cnv,as.numeric)
#colnames(cnv)<-cnv[1,]
#cnv<-cnv[-1,]














#a<-read.table("BIG.matrix2.txt")
#b<-read.csv("BIG.matrix2.csv")

#inter.genes<-intersect(TPM.mut.gene.intersect,colnames(cnv))
#inter.samples<-intersect(TPM.mut.sample.intersect,rownames(cnv))#

#BIG.matrix<-cbind(TPM.matrix[inter.samples,inter.genes],cnv[inter.samples,inter.genes],mut.aut[inter.samples,inter.genes])
#n.genes<-length(inter.genes)
#colnames(BIG.matrix)[1:n.genes]<-paste0("exp_",inter.genes)
#colnames(BIG.matrix)[(n.genes+1):(n.genes*2)]<-paste0("cnv_",inter.genes)
#colnames(BIG.matrix)[(n.genes*2+1):ncol(BIG.matrix)]<-paste0("mut_",inter.genes)
#g<-cbind(TPM.matrix[inter.samples,inter.genes],cnv[inter.samples,inter.genes])
# Read clinical file
##clinical<-read.delim("./clinical/clinical/Biotab/nationwidechildrens.org_clinical_patient_luad.txt",stringsAsFactors=FALSE)
#clinical<-clinical[3:nrow(clinical),c("bcr_patient_barcode","death_days_to")]
#rownames(clinical)<-NULL
#clinical<-arrange(clinical,bcr_patient_barcode)
#clinical
#Add clinical data to TPM.matrix

#x<-match(strtrim(rownames(BIG.matrix),12),clinical[,1])
#death_days_to<-clinical$death_days_to[x]
#BIG.matrix<-cbind(death_days_to,BIG.matrix)


#write.csv(BIG.matrix,"BIG.matrix.csv")
#write.csv(t(BIG.matrix),"tBIG.matrix.csv")

  
#Validate matrix
#rand.patient<-LUAD.files[sample(1:10,5, replace=F)]
#rand.gene<-sample(1:20531,5,replace=F)
#TPM.matrix[rand.gene,rand.patient]
#gene_id<-read.table(LUAD.files[1],header=T,stringsAsFactors=F)[,"gene_id"]
#split.gene.name<-t(as.data.frame(strsplit(gene_id[1:10],"|",fixed=T)))
#rownames(split.gene.name)<-NULL






#library("biomaRt")
#entrez=c("673","837")
# ensembl=useMart("ensembl")
#filters = listFilters(ensembl)
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#attributes = listAttributes(ensembl)
# goids=getBM(attributes=c('entrezgene','go_id'),values="672", mart=ensembl,fiters=filters)
 
#head(goids)


##big<-read.csv("BIG.matrix11.csv")
#colnames(big)[3:17926]<-paste0("exp_",inter.genes)
#write.csv(BIG.matrix,"BIG.matrix.csv")
#?write



#Adding missing samples from TPM with zeroes


#colnames(mut.cur)<-gsub(".","-",colnames(mut.cur),fixed=TRUE)
#Locations of genes in mut.cur that overlap gene_id
#gene.overlap<-match(gene_id[,"Symbol"],mut.cur[,1])
#sample.overlap<-match(strtrim(index$PatientID,15),colnames(mut.cur))
#Big.matrix<-cbind(TPM.matrix,t(mut.cur[which(gene.overlap!="NA"),which(sample.overlap!="NA")]))
