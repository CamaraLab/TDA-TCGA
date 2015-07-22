require(data.table)
anno_agilent<-read.delim("Google Drive/LAB on Drive/GBM/Agilent/agilent_annotations.txt",stringsAsFactors = F)
anno_agilent<-select(anno,EntrezGeneID,GeneSymbol) %>% distinct(EntrezGeneID)
index<-read.delim("Google Drive/LAB on Drive/GBM/Agilent/file_manifest.txt",stringsAsFactors = F,header = T)

files<-list.files("Google Drive/LAB on Drive/GBM/Agilent/Expression/",full.names = T)
matrix<-fread("Google Drive/LAB on Drive/GBM/Agilent/Expression/US14702406_251584710166_S01_GE2-v5_91_0806.txt_lmean.out.logratio.gene.tcga_level3.data.txt",stringsAsFactors = F,header=FALSE,data.table=F)  

samples<-matrix[1,2]
for (i in 2:length(files)) {
  data<-read.delim(files[i],stringsAsFactors = F,header=F)
  samples<-c(samples,data[1,2])
  matrix<-merge(matrix,data,by = "V1")
}

colnames(matrix)<-samples
rownames(matrix)<-matrix[,1]
matrix<-matrix[,-1]
matrix<-matrix[sort(rownames(matrix)),sort(colnames(matrix))]


rownames(matrix)<-paste0(rownames(matrix),"|",anno$EntrezGeneID[match(rownames(matrix),anno_agilent$GeneSymbol)])
rownames(matrix)[grepl("NA",rownames(matrix))]
matrix<-round(matrix,5)
rownames(matrix)
getwd()
write.csv(matrix,"Google Drive/LAB on Drive/GBM/Agilent/Agilent.Matrix.csv")
  
