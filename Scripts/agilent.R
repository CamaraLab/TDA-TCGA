require(data.table)
anno_agilent<-read.delim("Google Drive/LAB on Drive/GBM/Agilent/agilent_annotations.txt",stringsAsFactors = F)
anno_agilent<-select(anno,EntrezGeneID,GeneSymbol) %>% distinct(EntrezGeneID)
files<-list.files("Google Drive/LAB on Drive/GBM/Agilent/Expression/",full.names = T)
matrix<-fread("Google Drive/LAB on Drive/GBM/Agilent/Expression/US14702406_251584710166_S01_GE2-v5_91_0806.txt_lmean.out.logratio.gene.tcga_level3.data.txt",stringsAsFactors = F,header=T)  
samples<-colnames(matrix)[2]

for (i in 2:length(files)) {
  data<-fread(files[i],stringsAsFactors = F,header=T)
  samples<-c(samples,colnames(data)[2])
  matrix<-merge(matrix,data,by="Hybridization REF")
}




matrix<-as.matrix(matrix)
genes<-matrix[,1]
suppressWarnings(matrix<-apply(matrix,2,function(x) round (as.numeric(x),5)))
rownames(matrix)<-genes
matrix<-matrix[,-1]
matrix<-matrix[sort(rownames(matrix)),sort(colnames(matrix))]
rownames(matrix)<-paste0(rownames(matrix),"|",anno$EntrezGeneID[match(rownames(matrix),anno_agilent$GeneSymbol)])
matrix<-t(matrix)
#View(matrix)


write.csv(matrix,"Google Drive/LAB on Drive/GBM/Agilent/Agilent.Matrix.csv")
  
