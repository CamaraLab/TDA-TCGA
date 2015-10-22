# Environemnt settings
#setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
#library(org.Hs.eg.db)

suppressWarnings({
  suppressMessages ({
    library(dplyr,quietly = T,warn.conflicts = FALSE)
    library(stringr,quietly = T,warn.conflicts = FALSE)
    library(parallel,quietly = T,warn.conflicts = FALSE)
    library(getopt,quietly = T,warn.conflicts = FALSE)
  })
  
})


spec = matrix(c(
  "dir", "d",1, "character",
  "index", "i",1, "character",
  "anno", "a",1, "character",
  "anno_old", "o",1, "character",
  "project", "p", 1, "character",
  "cores", "q", 1, "integer"
  
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode


if ( is.null(arg$anno ) ) {arg$anno= "C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv"}
if ( is.null(arg$anno_old ) ) {arg$anno_old= "C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/anno_old_new.csv"}
if ( is.null(arg$cores ) ) {arg$cores= 1}


PROJECT_NAME<-arg$project

anno<-read.csv(arg$anno,as.is=T)
rownames(anno)<-anno[,1]
anno_old_new<-read.csv(arg$anno_old,as.is=T,row.names=1)


#Raw file list.
rsem_files<-list.files(arg$dir,full.names=T)
rsem_files<-rsem_files[grep(x=rsem_files,"rsem.genes.results")]

#Creating index file
index.file<-read.delim(arg$index,as.is=T)
index.cols<-c("Extract.Name","Comment..TCGA.Barcode.")
index<-unique(index.file[,index.cols]) # Remove duplicates
colnames(index)<-c("File","PatientID")
index<-arrange(index,File)

#Checking if index file matches rsem files
print ("Checking if index file matches rsem files (length)")
if (length(rsem_files)==nrow(index)) {
  print ("Length matches proceeding")
} else stop ("Number of samples in index does not match rsem files population")

#Creating gene list and file (gene_id)
gene_id_raw<-read.table(rsem_files[1],header=T,stringsAsFactors = F)[,"gene_id"] #Reading gene names from first file
gene_id_split<-strsplit(gene_id_raw,"|",fixed=TRUE)
symbol<-sapply(gene_id_split,"[[",1)
id<-as.numeric(sapply(gene_id_split,"[[",2))


#Replacing old ids with new
ids_to_replace_indices<-id %in% anno_old_new$Old_ID
ids_to_replace<-id[ids_to_replace_indices]
new_id_indices<-match(ids_to_replace,anno_old_new$Old_ID)
id[ids_to_replace_indices]<-anno_old_new$New_ID[new_id_indices]

symbol<-anno[as.character(id),"Symbol"] #Updating symbols avvording to Id's
gene_id<-as.data.frame(cbind(id,symbol),stringsAsFactors=F)

#Bring back withdrawn entrez ids
withdrawn_genes_indices<-as.numeric(rownames(gene_id[is.na(gene_id$id),]))
original_id<-sapply(gene_id_split[withdrawn_genes_indices],"[[",2)
gene_id[withdrawn_genes_indices,"id"]<-original_id
gene_id[withdrawn_genes_indices,"symbol"]<-rep("withdrawn",length(original_id))


# Reading scaled expression level into scale.estimates matrix
print (paste("Extracting TPM values from rsem files"))
print (paste("Utilizing",arg$cores,"CPU's cores"))
print ("This might take a few minutes")
cl <- makeCluster(arg$cores)
clusterExport(cl=cl, varlist=c("rsem_files"))
scale.estimates<-parSapply(cl,rsem_files,function (x) return(read.table(x,header=T)[,3]))
stopCluster(cl)
scale.estimates<-t(scale.estimates)




#Creating TPM matrix - 
print ("Creating TPM_matrix") 
samples<-index$PatientID
TPM.matrix<-as.data.frame(round(log2(1+scale.estimates*10^6),4))
rownames(TPM.matrix)<-samples
colnames(TPM.matrix)<-paste0(gene_id$symbol,"|",gene_id$id)
TPM.matrix<-TPM.matrix[sort(rownames(TPM.matrix)),sort(colnames(TPM.matrix))]
TPM.matrix<-round(TPM.matrix,5)
write.csv(TPM.matrix,paste0(PROJECT_NAME,"_Full_TPM_matrix.csv"))
print ("Done")
