#setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/SC-TDA/UDI")


KLD<-function(P,Q){
  
  #Measures Kullback–Leibler divergence between P and Q
  
  KL<-sum(P*log2i(P/Q))
}

JSD<-function(P,Q){
  #Measures JSD divergence between P and Q. JSD  =  0.5*KLD(P||M)+0.5*KLD(Q||M). 
  M<-0.5*(P+Q)
  jsd<-0.5*KLD(P,M)+0.5*KLD(Q,M)
  return(jsd)
}

log2i<-function(p) 
  # Special function that returns 0  if log2 argument is 0
  ifelse(log2(p)=="-Inf",0,p)

JSD1<-function(P,Q)
  
  0.5*sum(P*log2i(P))+0.5*sum(Q*log2i(Q))-0.5*sum((P+Q)*log2i((P+Q)/2))  

#Generating pii_file
cl <- makeCluster(as.numeric(arg$cores))
varlist=c("p_connectivity","arg","p_value","permute","c_vector","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns")
clusterExport(cl=cl, varlist=varlist,envir=environment())



pii_val<-read.csv("LUAD_Neigh_45_3_TPM.matrix.csv_pii_values_x.csv",row.names=FALSE)
columns_of_interest<-match(rownames(combined),colnames(pii_val))
pii<--pii_val[,columns_of_interest]
pi_val<--pii
pii_values<-read.csv("~/Desktop/Merging file/LUAD_Neigh_45_3_TPM.matrix.csv_pii_values.csv",row.names=1)
combined<-read.csv("~/Desktop/Merging file/combined.csv",row.names=1)

pii_values<-t(pii_values)
pii_values<-pii_values[rownames(x3),]
rownames(pii_values)
dim(pii_values)






genes_of_interest<-rownames(results3)[1:100] 
#genes_of_interest<-rownames(results[results$q.value<.05,])[1:100] 
pii_values<-pii_values_table(nodes,genes_of_interest)  

dim(pii_values)

JSD(pii_values[,2],pii_values[,3])


matrix1<-read.csv("TPM.matrix.LUAD.csv",row.names = 1)
results1<-read.csv("LUAD_Neigh_45_3_TPM.matrix.csv_results_rolling.csv",row.names = 1)
rownames(results)
#results

d<-ncol(pi_val)

a<-matrix(0,d,d)
for (i in 2:(d-1)) 
  for (j in (i+1):d)
      a[i,j]<-JSD1(x[,i],x[,j])


d<-1000
mylist.names <- 1:d
mylist <- as.list(rep(0, length(mylist.names)))


a<-matrix(0,d,d)

for (i in 1:(d-1))
  for (j in (i+1):d)
      a[i,j]<-identical(mat[,i],mat[,j])


for (i in 1:nrow(a))
{
          y<-which(a[i,]==1)
          mylist[[i]]<-c(i,y)
          a[,y]<-0
} 

sum(sapply(mylist,function(x) length(x)-1))




which 
    
sapply(seq_along(pii_values),function (x) sapply (pii_values[x],JSD) pii_values)

final<-read.csv("~/Desktop/Merging file/LUAD_Neigh_45_3_TPM.matrix.csv_results_final.csv",row.names=1)
rolling<-read.csv("~/Desktop/Merging file/LUAD_Neigh_45_3_TPM.matrix.csv_results_rolling.csv",row.names=1)
x<-cbind(final,rolling[genes,3:5])
x1<-x[x$pi_fraction!=1,]
x2<-x[x$q.value!<.05,]
x3<-x[intersect(rownames(x1),rownames(x2)),]
write.csv(x3,"~/Desktop/Merging file/combined.csv",row.names=TRUE)
combined<-read.csv("~/Desktop/Merging file/combined.csv",row.names=1)
colnames(pii_val)==rownames(combined)

columns_of_interest<-match(rownames(combined),colnames(pii_val))
pii<--pii_val[,columns_of_interest]
pi_val<--pii
pii_values<-read.csv("~/Desktop/Merging file/LUAD_Neigh_45_3_TPM.matrix.csv_pii_values.csv",row.names=1)
pii_values<-t(pii_values)
pii_values<-pii_values[rownames(x3),]
rownames(pii_values)
dim(pii_values)


y<-apply(a,1,function(x) which(x==1))



genes<-rownames(final)
x<-cbind(final,rolling[genes,3:5])
