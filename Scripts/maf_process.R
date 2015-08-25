library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(getopt)

spec = matrix(c(
  "maf", "m",1, "character",
  "anno", "a",1, "character",
  "anno_old", "o",1, "character",
  "project", "p", 1, "character"
  
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode


if ( is.null(arg$anno ) ) {arg$anno= "C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv"}
if ( is.null(arg$anno_old ) ) {arg$anno_old= "C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/anno_old_new.csv"}





PROJECT_NAME<-arg$project


anno<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv",as.is=T)
row.names(anno)<-anno[,1]
anno_old_new<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/anno_old_new.csv",as.is=T,row.names=1)




#wd<-paste0("c:/Users/Udi/Documents/TCGA-DATA/",PROJECT_NAME)
#setwd(wd)

#This file cleans and created a file wit


maf<-read.delim(arg$maf,header = TRUE,as.is=T,comment.char = "#")
colInfo<-c("Hugo_Symbol","Entrez_Gene_Id","Variant_Classification","Tumor_Sample_Barcode")
maf<-maf[,colInfo]

#Replacing unknown entrezids
#entrez_to_symbol<-as.list(org.Hs.egSYMBOL)
symbol_to_entrez<-as.list(org.Hs.egALIAS2EG)
symbol_to_entrez_short<-unlist(symbol_to_entrez[sapply(symbol_to_entrez,length)==1]) #Keeping only records with 1 entrez id

r<-which(maf$Entrez_Gene_Id==0) #0 indicates unknown in original maf
#View(maf[r,])
#length(r)
maf$Entrez_Gene_Id[r]<-symbol_to_entrez_short[maf$Hugo_Symbol[r]]
#View(maf[r,])
#sum(is.na(maf$Entrez_Gene_Id))
maf$Entrez_Gene_Id<-as.numeric(maf$Entrez_Gene_Id)
maf<-arrange(maf,Entrez_Gene_Id)
            


#Replacing old EntrezID's with new ones
print("Replacing old with new EntrezID's might take a few minutes")
ids_to_replace_indices<-maf$Entrez_Gene_Id %in% anno_old_new$Old_ID
ids_to_replace<-maf$Entrez_Gene_Id[ids_to_replace_indices]
new_id_indices<-match(ids_to_replace,anno_old_new$Old_ID)
maf$Entrez_Gene_Id[ids_to_replace_indices]<-anno_old_new$New_ID[new_id_indices]



#Sanity check - Check if unique EntrezID=Unique hugo_symbols
maf$Hugo_Symbol<-anno[as.character(maf$Entrez_Gene_Id),"Symbol"]
maf<-maf[complete.cases(maf),] #Removing unknown EntrezIDs/Hugo_symbols
print ("Unique id Sanity check:")
length(unique(maf$Entrez_Gene_Id))==length(unique(maf$Hugo_Symbol))

maf$Column_name<-paste0(maf$Hugo_Symbol,"|",maf$Entrez_Gene_Id) #Adding column name Symbol|Entrez for downstream traceback
#maf$Tumor_Sample_Barcode<-substring(maf$Tumor_Sample_Barcode,1,15) #Trim Fix sample name
maf$Synonymous<-(maf$Variant_Classification=="Silent")
maf$Exons_Length<-anno[match(maf$Entrez_Gene_Id,anno$EntrezID),"length"]  #Some length are unknown


#Creating synonymous and non-synonymous matrices
#Creating an empty matrices

all_genes<-sort(unique(maf$Column_name))
all_samples<-sort(unique(maf$Tumor_Sample_Barcode))
mat_syn<-matrix(0,length(all_samples),length(all_genes))
dimnames(mat_syn)<-list(all_samples,all_genes)
mat_non_syn<-mat_syn #Replicating mat_syn

#Creating table of Synonymous mutations for Sample vs Entrez_Gene_Id
t_syn<-with(maf[maf$Synonymous,],table(Tumor_Sample_Barcode,Column_name))
t_non_syn<-with(maf[!maf$Synonymous,],table(Tumor_Sample_Barcode,Column_name))

#Plugging tables into 0 matrices
mat_syn[rownames(t_syn),colnames(t_syn)]<-t_syn
mat_non_syn[rownames(t_non_syn),colnames(t_non_syn)]<-t_non_syn
mat_syn<-mat_syn[sort(rownames(mat_syn)),sort(colnames(mat_syn))]
mat_non_syn<-mat_non_syn[sort(rownames(mat_non_syn)),sort(colnames(mat_non_syn))]

#Binary matrix for connectivity score
mat_non_syn_bin<-ifelse(mat_non_syn>0,1,0) #Non synonymous binary matrix - will be used as input for c_score


#write.csv(mat_non_syn_bin,paste0("./Mutations/",PROJECT_NAME,"_Full_Mutations_binary.csv"))
#write.csv(mat_non_syn,paste0("./Mutations/",PROJECT_NAME,"_Full_Mutations_non_synonymous.csv"))
#write.csv(mat_syn,paste0("./Mutations/",PROJECT_NAME,"_Full_Mutations_synonymous.csv"))

print ("Writing mutation matrix files")
write.csv(mat_non_syn_bin,paste0(PROJECT_NAME,"_Full_Mutations_binary.csv"))
write.csv(mat_non_syn,paste0(PROJECT_NAME,"_Full_Mutations_non_synonymous.csv"))
write.csv(mat_syn,paste0(PROJECT_NAME,"_Full_Mutations_synonymous.csv"))

