clean_samples<-function(matrix) {
  samples<-substring(rownames(matrix),1,15)
  samples_to_keep_1<-which(sapply(samples, function (x) substring(x,14) %in% c("01","10","11")))
  #Removing duplicated records
  duplicated<-samples[duplicated(samples)]
  samples_to_keep_2<-which(!samples %in% duplicated)
  samples_to_keep<-samples[intersect(samples_to_keep_1,samples_to_keep_2)]
  matrix<-matrix[samples_to_keep,]
  rownames(matrix)<-samples_to_keep
  matrix<-matrix[sort(rownames(matrix)),sort(colnames(matrix))]
  return(as.matrix(matrix))
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
  #genes<-colnames(mat.A)
  sample.intersect<-intersect(samples,rownames(mat.B))
  #gene.intersect<-intersect(genes,colnames(mat.B))
  #mat<-mat.B[sample.intersect,gene.intersect]
  mat<-mat.B[sample.intersect,]
  
}

delta.zero.matrix<-function (mat.A,mat.B)
  # Creates delta matrix with missing rows from mut.A
  #nrow(mut.B+Delta)=mut.A - filled with zeros
{
  samples<-rownames(mat.A)
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




