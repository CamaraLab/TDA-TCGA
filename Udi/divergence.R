#setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/SC-TDA/UDI")


KLD<-function(P,Q){
  
  #Measures Kullback–Leibler divergence between P and Q
  LR<-log2(P/Q)
  KL<-sum(P*LR)
  
  return(KL)
}

JSD<-function(P,Q){
  #Measures JSD divergence between P and Q. JSD  =  0.5*KLD(P||M)+0.5*KLD(Q||M). 
  M<-0.5*(P+Q)
  jsd<-0.5*KLD(P,M)+0.5*KLD(Q,M)
  return(jsd)
}

JSD1<-function(P,Q)
  
  -0.5*sum(p*log2p)+
  JSDi<-function(p,q) -0.5*(p+q)*log2((p+q)/2)+0.5*p*log2p 
  jsd<-sapply(seq_along(P),)

#Generating pii_file

pii_values <-parSapply(cl,seq_along(columns),function (x) pii_calc (nodes,columns[x]))
pii_values<-cbind(1:length(nodes),pii_values)
colnames(pii_values)<-c("Node",colnames(matrix1)[columns])
print("Writing pii_values file:")
write.table(pii_values,paste0(arg$name,"_",arg$matrix,"_pii_values.csv"),sep=",",row.names=FALSE)


dim(pii_values)

JSD(pii_values[,2],pii_values[,3])


#matrix1<-read.csv("TPM.matrix.LUAD.csv",row.names = 1)
#results<-read.csv("LUAD_Neigh_45_3_TPM.matrix.csv_results_final.csv",row.names = 1)
#results

