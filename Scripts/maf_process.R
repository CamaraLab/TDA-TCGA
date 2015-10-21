suppressWarnings({
  suppressMessages ({
    library(data.table,quietly = T,warn.conflicts = FALSE)
    library(org.Hs.eg.db,quietly = T,warn.conflicts = FALSE)
    library(dplyr,quietly = T,warn.conflicts = FALSE)
    library(getopt,quietly = T,warn.conflicts = FALSE)
  })
  
})


spec = matrix(c(
  "maf", "m",1, "character",
  "anno", "a",1, "character",
  "anno_old", "o",1, "character",
  "project", "p", 1, "character"
  
), byrow=TRUE, ncol=4)


arg<-getopt(spec) #Conmment this line for debug mode

#arg<-list(project="BRCA",anno="Annotations.csv",anno_old_new="anno_old_new.csv",maf="BRCA_genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2_OCT-20-2015.maf")

if ( is.null(arg$anno ) ) {arg$anno= "C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv"}
if ( is.null(arg$anno_old ) ) {arg$anno_old= "C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/anno_old_new.csv"}



PROJECT_NAME<-arg$project
anno<-read.csv(arg$anno,as.is=T,row.names = 1)
anno$EntrezID<-rownames(anno)
anno_old_new<-read.csv(arg$anno_old,as.is=T,row.names=1)

#This file cleans and created a file wit
print ("Reading and reformatting maf file")
maf_file<-arg$maf
maf<-read.delim(maf_file,header = TRUE,as.is=T,comment.char = "#")
colInfo<-c("Hugo_Symbol","Entrez_Gene_Id","Variant_Classification","Tumor_Sample_Barcode")
maf<-maf[,colInfo]

#Replacing unknown entrezids (0) in original maf file
#entrez_to_symbol<-as.list(org.Hs.egSYMBOL)
print ("Annotations fix")
symbol_to_entrez<-as.list(org.Hs.egALIAS2EG)
symbol_to_entrez_short<-unlist(symbol_to_entrez[sapply(symbol_to_entrez,length)==1]) #Keeping only records with 1 entrez id
r<-which(maf$Entrez_Gene_Id==0) #0 indicates unknown in original maf

maf$Entrez_Gene_Id[r]<-symbol_to_entrez_short[maf$Hugo_Symbol[r]]
maf$Entrez_Gene_Id<-as.numeric(maf$Entrez_Gene_Id)
maf<-arrange(maf,Entrez_Gene_Id)




#Replacing old EntrezID's with new ones
print("Replacing old with new EntrezID's might take a few minutes")
ids_to_replace_indices<-maf$Entrez_Gene_Id %in% anno_old_new$Old_ID
ids_to_replace<-maf$Entrez_Gene_Id[ids_to_replace_indices]
new_id_indices<-match(ids_to_replace,anno_old_new$Old_ID)
maf$Entrez_Gene_Id[ids_to_replace_indices]<-anno_old_new$New_ID[new_id_indices]

#Replacing Hugo symbols with updated ones
removed_genes_for_log<-maf[!complete.cases(maf),]$Hugo_Symbol
maf$Hugo_Symbol<-anno[as.character(maf$Entrez_Gene_Id),"Symbol"]

#Removing unknown EntrezIDs/Hugo_symbols
maf<-maf[complete.cases(maf),] 

#Sanity check - Check if unique EntrezID=Unique hugo_symbols
print ("Unique id Sanity check:")
check<-length(unique(maf$Entrez_Gene_Id))==length(unique(maf$Hugo_Symbol))
print(check)
if (check==FALSE) {
  print ("Looking for ambivalentic ID/Symbol")
  unique_symbols<-unique(maf$Hugo_Symbol)
  x<-function(symbol) {
    y<-maf$Entrez_Gene_Id[maf$Hugo_Symbol==symbol]
    y<-length(unique(y))
  }
  z<-sapply(unique_symbols,x)
  genes_to_remove<-names(which(z>1))
  print ("Removing genes with duplicated ENTREZ Ids from the maf")
  print(genes_to_remove)
  removed_genes_for_log<-c(removed_genes_for_log,genes_to_remove)
  #Removing genes with duplicated Entrez ID
  rows_to_remove<-which(maf$Hugo_Symbol %in% genes_to_remove)
  maf<-maf[-rows_to_remove,]
  dim(maf)
  print ("Performing second sanity check:")
  check<-length(unique(maf$Entrez_Gene_Id))==length(unique(maf$Hugo_Symbol))
  print (check)
}



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

print ("Writing mutation matrix files")
write.csv(mat_non_syn_bin,paste0(PROJECT_NAME,"_Full_Mutations_binary.csv"))
write.csv(mat_non_syn,paste0(PROJECT_NAME,"_Full_Mutations_non_synonymous.csv"))
write.csv(mat_syn,paste0(PROJECT_NAME,"_Full_Mutations_synonymous.csv"))
write.table(maf,paste0("PROCESSED_",maf_file),sep="\t",row.names = FALSE)
cat("Genes removed in the process",file="outfile.txt",append=TRUE)
cat(removed_genes_for_log,file="outfile.txt",append=TRUE)