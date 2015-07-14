#This file cleans and created a file wit

maf<-read.table("../../../../../../Downloads/PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.automated.somatic.maf",header = TRUE,stringsAsFactors = FALSE)
head(maf)


#Replacing unknown entrezids
require(org.Hs.eg.db)
entrez_to_symbol<-as.list(org.Hs.egSYMBOL)
symbol_to_entrez<-as.list(org.Hs.egALIAS2EG)

symbol_to_entrez_short<-unlist(symbol_to_entrez[sapply(symbol_to_entrez,length)==1]) #Keeping only records with 1 entrez id
r<-which(maf$Entrez_Gene_Id==0) #0 indicates unknown in original maf
maf$Entrez_Gene_Id[r]<-symbol_to_entrez_short[maf$Hugo_Symbol[r]]

maf<-maf[order(maf$Entrez_Gene_Id),]
maf<-maf[complete.cases(maf),] #Removing unknown EntrezIDs

#Replacing old EntrezID's with new ones
old_new_id<-read.csv("../../../../../../Google Drive/Columbia/LAB/Rabadan/SC-TDA/Udi/entrez_oldid_conversion.csv")
for (i in 1:nrow(old_new_id)) 
  maf$Entrez_Gene_Id[maf$Entrez_Gene_Id==old_new_id[i,1]]<-old_new_id[i,2]


maf$Hugo_Symbol<-as.character(entrez_to_symbol[maf$Entrez_Gene_Id])
#Sanity check - Check if unique EntrezID=Unique hugo_symbols
length(unique(maf$Entrez_Gene_Id))==length(unique(maf$Hugo_Symbol))


maf$Tumor_Sample_Barcode<-substring(maf$Tumor_Sample_Barcode,1,15) #Trim Fix sample name
maf$Synonymous<-(maf$Variant_Classification=="Silent")
head(maf,30)

#Adding Exons_Length to file
anno<-read.csv("Annotations.csv")
maf$Exons_Length<-anno[match(maf$Entrez_Gene_Id,anno$EntrezID),"length"]  #Some length are unknown

#Creating synonymous and non-synonymous matrices
#Creating an empty matrices
mat_syn<-matrix(0,length(unique(maf$Tumor_Sample_Barcode)),length(unique(maf$Hugo_Symbol)))
rownames(mat_syn)<-unique(maf$Tumor_Sample_Barcode)
all_genes<-unique(maf$Hugo_Symbol)
colnames(mat_syn)<-all_genes
mat_non_syn<-mat_syn #Replicating mat_syn

#Creating table of Synonymous mutations for Sample vs Entrez_Gene_Id
t_syn<-with(maf[maf$Synonymous,],table(Tumor_Sample_Barcode,Hugo_Symbol))
t_non_syn<-with(maf[!maf$Synonymous,],table(Tumor_Sample_Barcode,Hugo_Symbol))

#Plugging tables into 0 matrices
mat_syn[rownames(t_syn),colnames(t_syn)]<-t_syn
mat_non_syn[rownames(t_non_syn),colnames(t_non_syn)]<-t_non_syn
mat_syn<-mat_syn[sort(rownames(mat_syn)),sort(colnames(mat_syn))]
mat_non_syn<-mat_non_syn[sort(rownames(mat_non_syn)),sort(colnames(mat_non_syn))]

#Output matrix for connectivity score
mat_non_syn_bin<-ifelse(mat_non_syn>0,1,0) #Non synonymous binary matrix - will be used as input for c_score


#Generating vector of genes lengths
exon_length_table<-cbind(maf$Hugo_Symbol,maf$Entrez_Gene_Id,maf$Exons_Length)
exon_length_table<-unique(exon_length_table)
rownames(exon_length_table)<-exon_length_table[,1]

Lg<-as.numeric(anno$length[match(colnames(mat_non_syn_bin),anno$Symbol)])
names(Lg)<-anno$Symbol[match(colnames(mat_non_syn_bin),anno$Symbol)]
genes_with_known_length<-names(Lg[!is.na(Lg)])
Lg<-Lg[genes_with_known_length] #Removing unknown length EntrezId's
L<-sum(Lg) #Total Coding region length
ns<-rowSums(mat_non_syn[,genes_with_known_length]) #Sum of non syn mutations for each row
#lambda_i<-Lg*ns/L


#G_scores type 1 - Based on  non-syn/(syn+non-syn) ratio
syn_ratio<-colSums(mat_non_syn)/(colSums(mat_syn)+colSums(mat_non_syn)) #G_Score
syn_ratio[syn_ratio=="NaN"]<-0 #Fixing 0 mutations columns 3

# G_Scores type 2 - Based on gene lengths
S_lambda<-rowSums(sapply(rownames(mat_non_syn),function (x) {
    lambda<-Lg*ns[x]/L 
    ans<-(-1)*mat_non_syn_bin[x,genes_with_known_length]*log(1-exp(-lambda))
}))

#G_Score 3
#Divides each element by the sum of the corresponding row sum.
# Returns zero in case of division by zero
  NS<-rowSums(mat_non_syn_bin)
  mat_NS<-mat_non_syn_bin/NS
  mat_NS[mat_NS=="NaN"]<-0
  NS<-colSums(mat_NS)

g_scores<-cbind(syn_ratio,S_lambda,NS)
write.csv(mat_non_syn_bin,"LUAD_Mutations_Binary.csv")
write.csv(mat_non_syn,"LUAD_Mutations_NS.csv")
write.csv(mat_syn,"LUAD_Mutations_S.csv")

require(rhdf5)
head(g_scores)
h5createFile("LUAD.h5")
H5close()
suppressWarnings({
  h5write(mat_non_syn_bin, "LUAD.h5","Mutations_Binary")
  h5write(mat_non_syn, "LUAD.h5","Mutations_NS")
  h5write(mat_syn, "LUAD.h5","Mutations_S")
  h5write(rownames(mat_non_syn_bin),"LUAD.h5","Mutations_Samples")
  h5write(colnames(mat_non_syn_bin),"LUAD.h5","Mutations_Genes")
})  
H5close()
h5ls("LUAD.h5")
