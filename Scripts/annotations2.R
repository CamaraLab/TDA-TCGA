history<-fread("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/gene_history_07_17_2015",data.table=FALSE,stringsAsFactors = F,na.strings = "-")

#Historical names handling
history<-history[history$V1=="9606",c(2,3,4)]
colnames(history)<-c("New_ID","Old_ID","Old_Name")

#Annotations file constrution
anno<-read.delim("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Homo_sapiens.gene_info_07_17_2015",skip=1,header=F,na.strings = "-",as.is=T)
anno<-anno[anno$V1==9606,c(2,3,5)]
colnames(anno)<-c("EntrezID","Symbol","Synonyms")
rownames(anno)<-anno[,1]

#Adding old ID reference from history file
anno$Old_ID<-history$Old_ID[match(anno$EntrezID,history$New_ID)]
anno$Old_Name<-history$Old_Name[match(anno$EntrezID,history$New_ID)]
symbol_conflicts<-filter(anno,duplicated(Symbol))$Symbol

head(anno)



#


h1<-history[complete.cases(history),]
h1[duplicated(h1$Old_ID_ID),]
head(history[(duplicated(history$New_ID[complete.cases(history)])),])
filter(h1,New_ID==5002)
head(anno)
head(history)
dim(anno)


cols<-c("refGene","Starts","Ends","EntrezID")
exons<-read.table("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/refseq_table4.txt",stringsAsFactors = F,header=T,comment.char="",col.names = cols)


head(exons)

