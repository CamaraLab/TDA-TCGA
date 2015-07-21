library(data.table)
library(org.Hs.eg.db)
library(rhdf5)
library(dplyr)

anno<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/Annotations.csv",as.is=T)
row.names(anno)<-anno[,1]
anno_old_new<-read.csv("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Annotations/anno_old_new.csv",as.is=T,row.names=1)



PROJECT_NAME<-"GBM"
wd<-paste0("c:/Users/Udi/Documents/TCGA-DATA/",PROJECT_NAME)
setwd(wd)

#This file cleans and created a file wit


maf<-read.delim("./Mutations/ucsc.edu_GBM.IlluminaGA_DNASeq_automated.Level_2.1.1.0.somatic.maf",header = TRUE,as.is=T,skip = 1)
colInfo<-c("Hugo_Symbol","Entrez_Gene_Id","Variant_Classification","Tumor_Sample_Barcode")
maf<-maf[,colInfo]

#Replacing unknown entrezids
#entrez_to_symbol<-as.list(org.Hs.egSYMBOL)
symbol_to_entrez<-as.list(org.Hs.egALIAS2EG)
symbol_to_entrez_short<-unlist(symbol_to_entrez[sapply(symbol_to_entrez,length)==1]) #Keeping only records with 1 entrez id

r<-which(maf$Entrez_Gene_Id==0) #0 indicates unknown in original maf
maf$Entrez_Gene_Id[r]<-symbol_to_entrez_short[maf$Hugo_Symbol[r]]

maf$Entrez_Gene_Id<-as.numeric(maf$Entrez_Gene_Id)
maf<-arrange(maf,Entrez_Gene_Id)
            


#Replacing old EntrezID's with new ones
for (i in 1:nrow(anno_old_new)) 
  maf$Entrez_Gene_Id[maf$Entrez_Gene_Id==anno_old_new$Old_ID[i]]<-anno_old_new$New_ID[i]



#Sanity check - Check if unique EntrezID=Unique hugo_symbols
maf$Hugo_Symbol<-anno[as.character(maf$Entrez_Gene_Id),"Symbol"]
maf<-maf[complete.cases(maf),] #Removing unknown EntrezIDs/Hugo_symbols
length(unique(maf$Entrez_Gene_Id))==length(unique(maf$Hugo_Symbol))

maf$Column_name<-paste0(maf$Hugo_Symbol,"|",maf$Entrez_Gene_Id) #Adding column name Symbol|Entrez for downstream traceback
maf$Tumor_Sample_Barcode<-substring(maf$Tumor_Sample_Barcode,1,15) #Trim Fix sample name
maf$Synonymous<-(maf$Variant_Classification=="Silent")

#Adding Exons_Length to file
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

#Output matrix for connectivity score
mat_non_syn_bin<-ifelse(mat_non_syn>0,1,0) #Non synonymous binary matrix - will be used as input for c_score

mat_non_syn_bin<-clean_samples(mat_non_syn_bin)
mat_non_syn<-clean_samples(mat_non_syn)
mat_syn<-clean_samples(mat_syn)


#Intersecting rows from TPM_matrix
TPM_matrix<-fread("C:/Users/Udi/Documents/TCGA-DATA/GBM/Expression/GBM_TPM_matrix.csv",data.table=F,select = 1)
samples_of_interest<-TPM_matrix[,1]
samples_of_interest<-intersect(samples_of_interest,rownames(mat_non_syn_bin))

mat_non_syn_bin<-mat_non_syn_bin[samples_of_interest,]
mat_non_syn<-mat_non_syn[samples_of_interest,]
mat_syn<-mat_syn[samples_of_interest,]
#Generating vector of genes lengths
#lambda_i<-Lg*ns/L


#G_scores type 1 - Based on  non-syn/(syn+non-syn) ratio
#syn_ratio<-colSums(mat_non_syn)/(colSums(mat_syn)+colSums(mat_non_syn)) #G_Score
#syn_ratio[syn_ratio=="NaN"]<-0 #Fixing 0 mutations columns 3

# G_Scores type 2 - Based on gene lengths
#S_lambda<-rowSums(sapply(rownames(mat_non_syn),function (x) {
#    lambda<-Lg*ns[x]/L 
#    ans<-(-1)*mat_non_syn_bin[x,genes_with_known_length]*log(1-exp(-lambda))
#}))
#supressWarnings(S_lambda<-as.numeric(c(S_lambda,setdiff(colnames(mat_non_syn),names(S_lambda)))))


#G_Score 3
#Divides each element by the sum of the corresponding row sum.
# Returns zero in case of division by zero
  #NS<-rowSums(mat_non_syn_bin)
  #mat_NS<-mat_non_syn_bin/NS
  #mat_NS[mat_NS=="NaN"]<-0
  #NS<-colSums(mat_NS)

#g_scores<-cbind(syn_ratio,S_lambda,NS)



write.csv(mat_non_syn_bin,paste0("./Mutations/",PROJECT_NAME,"_Mutations_binary.csv"))
write.csv(mat_non_syn,paste0("./Mutations/",PROJECT_NAME,"_Mutations_non_synonymous.csv"))
write.csv(mat_syn,paste0("./Mutations/",PROJECT_NAME,"_Mutations_synonymous.csv"))


h5file<-paste0("Mutations/",PROJECT_NAME,".h5")
h5createFile(h5file)
H5close()
suppressWarnings({
  h5write(mat_non_syn_bin, h5file,"_Mutations_Binary")
  h5write(mat_non_syn,h5file,"Mutations_NS")
  h5write(mat_syn,h5file,"Mutations_S")
  h5write(rownames(mat_non_syn_bin),h5file,"Mutations_Samples")
  h5write(colnames(mat_non_syn_bin),h5file,"Mutations_Genes")
})  

print(h5ls(h5file))
H5close()
