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
varlist=c("p_connectivity","arg","p_value","permute","c_vector","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns","genes_of_interest")
clusterExport(cl=cl, varlist=varlist,envir=environment())



genes_of_interest<-rownames(results[results$q.value<.05,])  
pii_values<-pii_values_table(nodes,genes_of_interest)  

dim(pii_values)

JSD(pii_values[,2],pii_values[,3])


matrix1<-read.csv("TPM.matrix.LUAD.csv",row.names = 1)
results<-read.csv("LUAD_Neigh_45_3_TPM.matrix.csv_results_final.csv",row.names = 1)
rownames(results)
#results

