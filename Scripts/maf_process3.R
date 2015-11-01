suppressWarnings({
  suppressMessages ({
    library(data.table,quietly = T,warn.conflicts = FALSE)
    library(dplyr,quietly = T,warn.conflicts = FALSE)
    library(getopt,quietly = T,warn.conflicts = FALSE)
    library(parallel,quietly = T,warn.conflicts = FALSE)
  })
  
})


spec = matrix(c(
  "tar", "f",1, "character",
  "project", "p", 1, "character"
  
), byrow=TRUE, ncol=4)


arg<-getopt(spec) #Conmment this line for debug mode

arg<-list(project="BLCA",tar="gdac.broadinstitute.org_BLCA.Mutation_Packager_Oncotated_Raw_Calls.Level_3.2015082100.1.0.tar.gz")



PROJECT_NAME<-arg$project

#This file cleans and created a file wit


print ("Extracting oncotator files:")
extraction_dir<-"Extracted_Onconator"
dir.create(extraction_dir)

untar(tarfile = arg$tar,compressed = "gzip",exdir = extraction_dir)


print ("Reading oconator files")
mut_files_names<-list.files(extraction_dir,full.names = T,recursive = T,pattern = "TCGA")

colInfo<-c("Hugo_Symbol","Entrez_Gene_Id","Variant_Classification","Variant_Type","Tumor_Sample_Barcode")

print("Constructing maf file")
maf<-NULL
for (file in mut_files_names) {
  mut_file<-read.delim(file,header = TRUE,as.is=T,comment.char = "#")    
  maf<-rbind(maf,mut_file[,colInfo])
}

print ("Deleting extraction dir")
unlink(extraction_dir,recursive = T,force = T)

r<-which(maf$Entrez_Gene_Id==0) #0 indicates unknown in original mut file
removed_rows_for_log<-maf$Hugo_Symbol[maf$Entrez_Gene_Id==0]
removed_rows_for_log<-as.data.frame(table(removed_rows_for_log))
removed_rows_for_log<-arrange(removed_rows_for_log,desc(Freq))
removed_rows_for_log<-cbind(removed_rows_for_log,0)
colnames(removed_rows_for_log)<-c("Hugo_Symbol","Occurences","Entrez_Gene_Id")

maf<-maf[-r,]

maf<-arrange(maf,Entrez_Gene_Id)




#Replacing old EntrezID's with new ones
#print("Replacing old with new EntrezID's might take a few minutes")
#ids_to_replace_indices<-maf$Entrez_Gene_Id %in% anno_old_new$Old_ID
#ids_to_replace<-maf$Entrez_Gene_Id[ids_to_replace_indices]
#new_id_indices<-match(ids_to_replace,anno_old_new$Old_ID)
#maf$Entrez_Gene_Id[ids_to_replace_indices]<-anno_old_new$New_ID[new_id_indices]

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
  x1<-function(entrez) {
    y<-maf$Hugo_Symbol[maf$Entrez_Gene_Id==entrez]
    y<-length(unique(y))
  }
  
  cl <- makeCluster(4)
  clusterExport(cl=cl,varlist=c("maf"),envir=environment())
  
  
  z<-parSapply(cl,unique_symbols,x)
  z1<-parSapply(cl,unique(maf$Entrez_Gene_Id),x1)
  print ("non unique corresponding entrez id for genes:")
  genes_to_remove<-names(which(z>1))
  print(genes_to_remove)
  print ("non unique corresponding genes symbol  for entrez:")
  entrez_to_remove<-unique(maf$Entrez_Gene_Id)[z1>1]
  print(entrez_to_remove)
  conflicted_genes_entrez<-c(genes_to_remove,entrez_to_remove)
  #Removing genes with duplicated Entrez ID
  #rows_to_remove<-which(maf$Hugo_Symbol %in% genes_to_remove)
  #maf<-maf[-rows_to_remove,]
  #dim(maf)
  #print ("Performing second sanity check:")
  #check<-length(unique(maf$Entrez_Gene_Id))==length(unique(maf$Hugo_Symbol))
  #print (check)
}

#removed_genes_for_log<-(as.data.frame(sort(table(removed_genes_for_log),decreasing=TRUE)))
#colnames(removed_genes_for_log)<-"Events"

maf$Column_name<-paste0(maf$Hugo_Symbol,"|",maf$Entrez_Gene_Id) #Adding column name Symbol|Entrez for downstream traceback
#maf$Tumor_Sample_Barcode<-substring(maf$Tumor_Sample_Barcode,1,15) #Trim Fix sample name
maf$Synonymous<-(maf$Variant_Classification=="Silent")
#maf$Exons_Length<-anno[match(maf$Entrez_Gene_Id,anno$EntrezID),"length"]  #Some length are unknown


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
write.table(maf,paste0("PROCESSED_MAF_",PROJECT_NAME,"_",Sys.Date(),".maf"),sep="\t",row.names = FALSE)
write.csv(removed_rows_for_log,"removed_rows.csv")


