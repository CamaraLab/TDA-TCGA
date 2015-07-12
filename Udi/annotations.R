library(GenomicFeatures)
symbols<-c("STAT1","CXCL10","ACTB","PDCD1")
symbols<-"all"
gene_length<-hg19GeneLengths(symbols)
gene_length

#Parsing gene/exons/length table
exons<-read.delim("../Downloads/refseq_table.txt",stringsAsFactors = F,header=T)
exons$EntrezID<-as.numeric(unlist(strsplit(exons$EntrezID,",")))
exons<-subset(exons,!is.na(EntrezID)) #Removing undetected EntrezID's

exons$Starts<-sapply(exons$Starts,function (x) strsplit(x,","))
exons$Ends<-sapply(exons$Ends,function (x) strsplit(x,","))

entrezid<-unique(exons$EntrezID)
gene_length<-sapply(entrezid,function (id){
  starts<-as.numeric(unlist(exons$Starts[exons$EntrezID==id]))
  ends<-as.numeric(unlist(exons$Ends[exons$EntrezID==id]))
  range<-reduce(IRanges(starts,ends))
  ans<-sum(width(range))
})
names(gene_length)<-as.character(entrezid)
gene_length


GRangesList(1,t2)



maf<-read.table("../Downloads/PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.automated.somatic.maf",header = TRUE,stringsAsFactors = FALSE)


entrez_to_symbol<-as.list(org.Hs.egSYMBOL)
symbol_to_entrez<-as.list(org.Hs.egALIAS2EG)



#Replacing unknown entrezids
require(org.Hs.eg.db)
xx<-as.list(org.Hs.egALIAS2EG)
r<-which(maf$Entrez_Gene_Id==0)
for (i in r) 
  if (length(xx[[maf$Hugo_Symbol[i]]])==1)
  {maf$Entrez_Gene_Id[i]<-as.character(xx[maf$Hugo_Symbol[i]])}

maf<-maf[order(maf$Entrez_Gene_Id),]
head(maf)
maf<-subset(maf,Entrez_Gene_Id!=0)
head(maf)

#Replacing old EntrezID's with new ones
old_new_id<-read.csv("../Google Drive/Columbia/LAB/Rabadan/SC-TDA/Udi/entrez_oldid_conversion.csv")
for (i in 1:nrow(old_new_id)) 
maf$Entrez_Gene_Id[maf$Entrez_Gene_Id==old_new_id[i,1]]<-old_new_id[i,2]

maf$Tumor_Sample_Barcode<-substring(maf$Tumor_Sample_Barcode,1,15) #Trim Fix sample name
maf$Hugo_Symbol<-entrez_to_symbol[maf$Entrez_Gene_Id]
maf$Synonymous<-(maf$Variant_Classification=="Silent")
head(maf,30)
maf$Exons_Length<-gene_length[maf$Entrez_Gene_Id] 

samples<-unique(maf$Tumor_Sample_Barcode)
genes<-unique(maf$Hugo_Symbol)
#Check if unique EntrezID!=Unique hugo_symbols
length(unique(maf$Entrez_Gene_Id))==length(unique(maf$Hugo_Symbol))
ans<-sapply(symbol,function(x) length(unique(z[z$NAME==as.character(x),1])))

#Creating an empty matrices
mat_syn<-matrix(0,length(samples),length(unique(maf$Hugo_Symbol)))
rownames(mat_syn)<-samples
colnames(mat_syn)<-unique(maf$Entrez_Gene_Id)
mat_non_syn<-mat_syn

#Creating table of Synonymous mutations for Sample vs Entrez_Gene_Id
t_syn<-with(maf[maf$Synonymous,],table(Tumor_Sample_Barcode,Entrez_Gene_Id))
t_non_syn<-with(maf[!maf$Synonymous,],table(Tumor_Sample_Barcode,Entrez_Gene_Id))

#Plugging tables into matrices
mat_syn[rownames(t_syn),colnames(t_syn)]<-t_syn
mat_non_syn[rownames(t_non_syn),colnames(t_non_syn)]<-t_non_syn
colnames(mat_syn)<-entrez_to_symbol[colnames(mat_syn)]
colnames(mat_non_syn)<-entrez_to_symbol[colnames(mat_non_syn)]

mat_syn["TCGA-17-Z031-01",1:5]
