
log2i<-function(p)  {
  # Special function that returns 0  if log2 argument is 0
  w<-log2(p)
  ifelse(w=="-Inf",0,w)
}

#JSD1<-function (P,Q) {
#  return(sqrt(sum(0.5*P*log2i(P)+0.5*Q*log2i(Q)-0.5*(P+Q)*log2i(0.5*(P+Q)))))
#}

JSD2<-function(P,Q) {
  P1<-replace(P,P==0,1) #All Values of 1 will get log2 ==0
  Q1<-replace(Q,Q==0,1)
  return(sqrt(sum(0.5*P1*log2(P1)+0.5*Q1*log2(Q1)-0.5*(P+Q)*log2i(0.5*(P+Q))))) 
}

columns_cutoff<-function (final_results,q_value_cutoff=.05,min_frac=0,max_frac=1,equal=FALSE){
 #Recives final_results matrix,q_value,min_frac, and max_frac. And retreives columns of interest
 
  if (equal) {
    final_results<-subset(final_results,q_value<=q_value_cutoff)
    final_results<-subset(final_results,pi_fraction<=max_frac)
    final_results<-subset(final_results,pi_fraction>=min_frac)  
  } else{
    final_results<-subset(final_results,q_value<q_value_cutoff)
    final_results<-subset(final_results,pi_fraction<max_frac)
    final_results<-subset(final_results,pi_fraction>min_frac)
  }
  
  genes_of_interest<-rownames(final_results)
  
}

JSD_matrix<-function (pii_values,genes_of_interest,cores,d=0,fast=TRUE) {
  
  pii_values<-pii_values[,-1] #Subsetting pii matrix
  if (d==0) d<-length(genes_of_interest)
  
  if (!fast) {
    print ("Building matrix the slow way (no parallel):")
    mat<-matrix(0,d,d)
    for (i in 1:(d-1)) {
      print (paste("working on row",i,"of",d))
      for (j in (i+1):d) {
        mat[i,j]<-JSD2(pii_values[,i],pii_values[,j])
        mat[j,i]<-mat[i,j]
      }
    }
        
  } else {
    
    print ("Preparing parallel environment")
    
    cl <- makeCluster(cores)
    pool<-clusterSplit(cl,1:d)
    varlist=c("JSD2","d","pool","pii_values","log2i")
    clusterExport(cl=cl, varlist=varlist,envir=environment())
    
    print("Building matrix the fast way")
    write.table(0,"row.csv",col.names=FALSE,row.names=FALSE)
    ans<-parLapply(cl,seq_along(pool),function (x) {
    #ans<-lapply(seq_along(pool), function(x) {
      a1<-matrix(0,d,length(pool[[x]])) #Auxilary matrix split big matrix into cores size matrices
      for (j in min(pool[[x]]):max(pool[[x]])) {
        write.table(j,"row.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
        for (i in (1:d)) {
          
          if (i>=j) {
            cols<-j-(x-1)*length(pool[[1]])
            a1[i,cols]<-JSD2(pii_values[,i],pii_values[,j])
            #a1[cols,i]<-a1[i,cols]}
            
          }
        }
      }
        
      return(a1)
    })
    
    mat<-NULL
    for (i in 1:length(ans)) mat<-cbind(mat,ans[[i]])
    for (i in 1:(d-1))
      for (j in (i+1):d) mat[i,j]<-mat[j,i]
  
  }
  print("done")
  return(mat)
  
}
results<-read.csv("E:/Merging2/LUAD_Neigh_45_3_TPM.matrix.csv_results_final-3.csv",row.names=1)
genes_of_interest<-columns_cutoff(results,q_value_cutoff=.05,min_frac=-1,max_frac=1,equal=FALSE)
genes_of_interest

#pii1<-as.matrix(fread("E:/Merging2/LUAD_Neigh_45_3_TPM.matrix.csv-192553-2015-06-17_pii_values.csv",data.table=FALSE))
#pii1<-pii1[,-1]

e_slow<-JSD_matrix(pii_values = pii1,genes_of_interest = genes_of_interest,cores = 2,d = 25,fast =FALSE)
e_fast<-JSD_matrix(pii_values = pii1,genes_of_interest = genes_of_interest,cores = 2,d = 0,fast =TRUE)
  
dim(e_slow)
dim(e_fast)
identical(e_slow,e_fast)






  
  
  
  
