KLD<-function(P,Q){
  
  #Measures Kullback–Leibler divergence between P and Q
  P<-P/sum(P)
  Q<-Q/sum(Q)
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

pii_calc<-function(nodes,column)
{  #Calculate pi value(mean) for a particular column in a particular node.
  # Input: nodes list and one specific column(i) of interest, output: pi
  if (arg$log2==TRUE) {
    ei<-sapply(nodes,function (x) log2(1+mean(matrix1[x+1,column]))) 
    
  } else ei<-sapply(nodes,function (x) mean(matrix1[x+1,column]))
  
  total_ei<-sum(ei)  
  if (total_ei==0) {pii<-ei} else pii<-ei/total_ei #Preventing diviosin by zero
  return(pii)
}



matrix1<-read.csv("TPM.matrix.LUAD.csv",row.names = 1)
results<-read.csv("LUAD_Neigh_45_3_TPM.matrix.csv_results_final.csv",row.names = 1)
results




ans<-sapply(rownames(results)[1:4],function (x) {
  pii<-pii_calc(nodes,x)
  pii_mean<-mean(pii)
  pii_sd<-sd(pii)
  pii_frac<-sum(pii!=0)/length(pii)
  ans<-cbind(pii_mean,pii_sd,pii_frac)
  return(ans)
})
ans
t(ans)

