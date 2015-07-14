#This file will get table from ucsc table browser and extract coding region length
#Considering different transcripts for same gene by taking unioned length
#Output is a table with entrezId-gene,length
library(GenomicFeatures)
require(org.Hs.eg.db)
entrez_to_symbol<-as.list(org.Hs.egSYMBOL)

#Parsing gene/exons/length table
cols<-c("refGene","Starts","Ends","EntrezID")
exons<-read.table("../../../../../../Downloads/refseq_table4.txt",stringsAsFactors = F,header=T,comment.char="",col.names = cols)
#exons$EntrezID<-as.numeric(unlist(strsplit(exons$EntrezID,","))) #Cleaning EntrezID chracters
exons<-subset(exons,!is.na(EntrezID)) #Removing undetected EntrezID's
exons$Starts<-sapply(exons$Starts,function (x) strsplit(x,",")) #Extracting and cleaning starts positions
exons$Ends<-sapply(exons$Ends,function (x) strsplit(x,",")) #Extracting and cleaning ends positions

#Calculationg unioned length for each entrezID
entrezid<-unique(exons$EntrezID)
gene_length<-sapply(entrezid,function (id){
  starts<-as.numeric(unlist(exons$Starts[exons$EntrezID==id]))
  ends<-as.numeric(unlist(exons$Ends[exons$EntrezID==id]))
  range<-reduce(IRanges(starts,ends)) #Removing overlapping regions
  ans<-sum(width(range)) #Calculating total length
})
names(gene_length)<-as.character(entrezid) #Matching gene_length to calculated length
symbols<-as.character(entrez_to_symbol[names(gene_length)])
anno<-as.data.frame(cbind(names(gene_length),symbols,gene_length),stringsAsFactors = F)
colnames(anno)<-c("EntrezID","Symbol","length")
anno<-anno[order(as.numeric(anno$EntrezID)),]

write.csv(anno,"Annotations.csv",row.names=F)

