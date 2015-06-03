get_symbol<-function (unknown_id) 
  # Gets ENTREZID and return Symbol, if not exist return "?"
  {
  unknown_id<-as.character(unknown_id)
  allkeys<-keys(org.Hs.eg.db,keytype="ENTREZID")
  if (unknown_id %in% allkeys) 
    unknown_id<-select(org.Hs.eg.db,unknown_id,keytype="ENTREZID",columns="SYMBOL")$SYMBOL
    else unknown_id="?"
  }

fix_symbols<-function (gene_id_raw)
  # Gets gene_id with symbols and fixes the unknown ones using annotation library.
{
  gene_id_split<-strsplit(gene_id_raw,"|",fixed=TRUE)
  #rownames(gene_id_split)<-gene_id_split[,2]
  #gene_id<-sapply(gene_id_split,1,function (x) if (x[1]=="?") x[1]<-get_symbol(x[2]) else x[1]<-x[1])
  symbol<-sapply(gene_id_split,function (x) if (x[1]=="?") x[1]<-get_symbol(x[2]) else x[1]<-x[1])
  id<-sapply(gene_id_split,"[[",2)
  gene_table<-cbind(symbol,id)
  }

intersect.mat<-function(mat.A,mat.B)
  #mat.A is the base matrix, mat.A will be bigger
  #mat.B is intersection of genes and samples
{
  samples<-rownames(mat.A)
  genes<-colnames(mat.A)
  sample.intersect<-intersect(samples,rownames(mat.B))
  gene.intersect<-intersect(genes,colnames(mat.B))
  mat<-mat.B[sample.intersect,gene.intersect]
  
}

delta.zero.matrix<-function (mat.A,mat.B)
  # Creates delta matrix with missing rows from mut.A
  #nrow(mut.B+Delta)=mut.A - filled with zeros
{
  samples<-rownames(mat.A)
  genes<-colnames(mat.A)
  missing_samples<-setdiff(samples,rownames(mat.B))
  delta<-matrix(0,nrow(mat.A)-nrow(mat.B),ncol(mat.B))
  rownames(delta)<-missing_samples
  return(delta)
  }


fix_patient_id<-function(matrix1)
  #input-matrix,Change . to - in patinet_ID, patient-ID needs to be in the rows
{ rownames(matrix1)<-gsub(".","-",rownames(matrix1),fixed=TRUE)
  return(matrix1)}



remove.after<-function(string,character)
#Remove everything after special character from string
{
  if (grepl("|",string,fixed=TRUE)==TRUE){
    b<-unlist(strsplit(string,""))
    loc<-which(b==character)
    a<-b[1:(loc-1)]
    ans<-paste(a,collapse="")}
    
  if (grepl("|",string,fixed=TRUE)==FALSE) ans<-string
 }   







m_join<-function (matrix1,matrix2,rows,cols)
  matrix2<-sample.intersect<-intersect(rownames(TPM.matrix),colnames(mut.aut))
  gene.intersect<-intersect(gene_id_cur[,"Symbol"],mut.aut[,1])
  BIG.matrix<-cbind(TPM.matrix[sample.intersect,gene.intersect]
#Read mutation data
mut.aut<-read.delim("./MUTATIONS/UCSC/genomicMatrix.automated",stringsAsFactors=FALSE)
#mut.cur<-read.delim("./MUTATIONS/UCSC/genomicMatrix.curated",stringsAsFactors=FALSE)
mut.aut<-as.data.frame(mut.aut,stringsAsFactors=FALSE)
colnames(mut.aut)<-gsub(".","-",colnames(mut.aut),fixed=TRUE)
#colnames(mut.cur)<-gsub(".","-",colnames(mut.cur),fixed=TRUE)
#Locations of genes in mut.cur that overlap gene_id
#gene.overlap<-match(gene_id[,"Symbol"],mut.cur[,1])
#sample.overlap<-match(strtrim(index$PatientID,15),colnames(mut.cur))
#Big.matrix<-cbind(TPM.matrix,t(mut.cur[which(gene.overlap!="NA"),which(sample.overlap!="NA")]))

sample.intersect<-intersect(rownames(TPM.matrix),colnames(mut.aut))
gene.intersect<-intersect(gene_id_cur[,"Symbol"],mut.aut[,1])
BIG.matrix<-cbind(TPM.matrix[sample.intersect,gene.intersect],t(mut.aut[gene.intersect,sample.intersect]))
colnames(BIG.matrix)[17925:35848]<-gene.intersect





r<-fix_symbols(gene_id_split)
}

