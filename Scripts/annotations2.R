library(data.table)
library(GenomicFeatures)
setwd("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations")
#Historical names handling
history<-fread("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/gene_history_07_17_2015",data.table=FALSE,stringsAsFactors = F,na.strings = "-")
history<-history[history$V1=="9606",c(2,3,4)]
colnames(history)<-c("New_ID","Old_ID","Old_Name")
rownames(history)<-history$Old_ID
write.csv(history,"anno_old_new.csv")

#Gene length from refseq
cols<-c("refGene","Starts","Ends","EntrezID")
exons<-read.table("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/refseq_table4.txt",stringsAsFactors = F,header=T,comment.char="",col.names = cols)
#exons$EntrezID<-as.numeric(unlist(strsplit(exons$EntrezID,","))) #Cleaning EntrezID chracters
exons<-subset(exons,!is.na(EntrezID)) #Removing undetected EntrezID's
exons$Starts<-sapply(exons$Starts,function (x) strsplit(x,",")) #Extracting and cleaning starts positions
exons$Ends<-sapply(exons$Ends,function (x) strsplit(x,",")) #Extracting and cleaning ends positions

#Calculationg unioned length for each entrezID
print("Calculating gene lengths, it might take a few minutes")
entrezid<-unique(exons$EntrezID)
gene_length<-sapply(entrezid,function (id){
  starts<-as.numeric(unlist(exons$Starts[exons$EntrezID==id]))
  ends<-as.numeric(unlist(exons$Ends[exons$EntrezID==id]))
  range<-reduce(IRanges(starts,ends)) #Removing overlapping regions
  ans<-sum(width(range)) #Calculating total length
})
names(gene_length)<-as.character(entrezid) #Matching gene_length to calculated length


#Annotations file constrution
anno<-read.delim("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Homo_sapiens.gene_info_07_17_2015",skip=1,header=F,na.strings = "-",as.is=T)
anno<-anno[anno$V1==9606,c(2,3,5)]
colnames(anno)<-c("EntrezID","Symbol","Synonyms")
rownames(anno)<-anno[,1]


#symbol_conflicts<-filter(anno,duplicated(Symbol))$Symbol

#Adding Exonic length
anno$length<-gene_length[match(anno$EntrezID,names(gene_length))]


write.csv(anno,"Annotations2.csv",row.names=F)
