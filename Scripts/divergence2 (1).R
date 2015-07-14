library(parallel)
library(data.table)

log2i<-function(p)  {
  # Special function that returns 0  if log2 argument is 0
  w<-log2(p)
  ifelse(w=="-Inf",0,w)
}


JSD2<-function(P,Q) {
  #Calculates JSD distance between distribution P and Q
  P1<-replace(P,P==0,1) #All Values of 1 will get log2 ==0
  Q1<-replace(Q,Q==0,1)
  return(sqrt(sum(0.5*P1*log2(P1)+0.5*Q1*log2(Q1)-0.5*(P+Q)*log2i(0.5*(P+Q))))) 
}

columns_cutoff<-function (final_results,q_value_cutoff=.05,min_frac=0,max_frac=1,equal=FALSE){
 #Recives final_results matrix,q_value,min_frac, and max_frac. And retreives columns of interest
  parameters<<-paste("q_value_cut_off=",q_value_cutoff," ",min_frac, " <frac< ",max_frac, "equal=",equal)
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

JSD_matrix<-function (mat1,mat2,cores,d=0,fast=TRUE) {
  #This function takes  matrix and can do JSD calculation on it paralleli by splitting
  #The matrix to cores size lements and then rejoin. This function also has slow mode with no parallel
  
 genes_of_interest_1<-intersect(colnames(mat1),genes_of_interest_1)
 genes_of_interest_2<-intersect(colnames(mat2),genes_of_interest_2)
  mat1<-mat1[,genes_of_interest_1] # Subsetting for genes of interest
  mat2<-mat2[,genes_of_interest_2] # Subsetting for genes of interest
  
  
  ##Initializing log file
  log_file<-paste0("JSD_matrix-",Sys.Date(),".log")
  write.table(parameters,log_file,col.names=FALSE)
  write.table(results_file_location,log_file,append = TRUE,col.names=FALSE)
  write.table(pii_values_file_location,log_file,append = TRUE,col.names=FALSE)
  
  #pii_values<-pii_values[,-1] #Subsetting pii matrix
  if (d==0) {
    d1<-length(genes_of_interest_1)
    d2<-length(genes_of_interest_2)
  } else { d1<-d
          d2<-length(genes_of_interest_2)
          }
  if (!fast) {
    print ("Building matrix the slow way (no parallel):")
    mat<-matrix(0,d,d)
    #if (!identical(mat1[,genes_of_interest_1]==mat2[,genes_of_interest_2])) { #In this case we canc alculate only half matrix minus diagonal
      if (f!=1) {
      for (i in 1:(d-1)) {
        print (paste("working on row",i,"of",d1))
        write.table(paste0(Sys.time(),i),log_file,append=TRUE,col.names=FALSE,row.names=FALSE)
        for (j in (i+1):d) {
          mat[i,j]<-JSD2(mat1[,i],mat1[,j])
          mat[j,i]<-mat[i,j]
        }
      }  
    } else {  #In this case all columns need to calculated against all
      mat<-matrix(0,d1,d2)  
      for (i in 1:d1) {
        print (paste("working on col",i,"of",d1))
        write.table(paste0(Sys.time(),i),log_file,append=TRUE,col.names=FALSE,row.names=FALSE)
        for (j in 1:d2) {
          #print(i)
          #print(j)
          mat[i,j]<-JSD2(mat1[,i],mat2[,j])
          #mat[j,i]<-mat[i,j]
        }
      }
    
        
  }
    } else {
    
    print ("Preparing parallel environment")
    
    
    l<-list() # Preapring list that will handle cores jobs
    l<-as.list(rep(NA,cores) )
    
    element_size<-ceiling(d/cores)
    for (i in 1:cores)
    { 
      for (j in 1:element_size)
        l[[i]][j]<-(1:d)[(i-1)*element_size+j] }
    l[[cores]]<-l[[cores]][!is.na(l[[cores]])]     
    pool<-l
  
    cl <- makeCluster(cores)  
    varlist=c("log_file","JSD2","d","d1","d2","pool","log2i","mat1","mat2","f")
    clusterExport(cl=cl, varlist=varlist,envir=environment())
    
    print("Building matrix the fast way")
    
    ans<-parLapply(cl,seq_along(pool),function (x) {
    #ans<-lapply(seq_along(pool),function (x) {
      a1<-matrix(0,d2,length(pool[[x]])) #Auxilary matrix split big matrix into cores size matrices
      for (j in min(pool[[x]]):max(pool[[x]])) {
        write.table(paste0(Sys.time(),j),log_file,append=TRUE,col.names=FALSE,row.names=FALSE)
        for (i in (1:d2)) {
          if (f!=1) {
            if (i>=j) {
              cols<-j-(x-1)*length(pool[[1]])
              a1[i,cols]<-JSD2(mat1[,i],mat2[,j])
              
            } 
          } else {
            cols<-j-(x-1)*length(pool[[1]])
            a1[i,cols]<-JSD2(mat1[,j],mat2[,i])
          
            
          }
        }
      }
        
      return(a1)
    })
    print("Releasing CPU cores")
    stopCluster(cl)
    
    mat<-NULL #Combining splited matrices and mirroring values
    if (f!=1) {
      
      for (i in 1:(d-1))
        for (j in (i+1):d) mat[i,j]<-mat[j,i]  
    } else for (i in 1:length(ans)) mat<-cbind(mat,ans[[i]])
    
  
  }
  
  print("done")
  rownames(mat)<-genes_of_interest_2[1:nrow(mat)]
  colnames(mat)<-genes_of_interest_1[1:ncol(mat)]
  return(mat)
  
}

results_file_location<<-"LUAD_Neigh_45_3_Mut.matrix.csv-180908-2015-06-17_results_final.csv"
pii_values_file_location<<-"LUAD_Neigh_45_3_Mut.matrix.csv-180908-2015-06-17_pii_values.csv"

results1<-read.csv(results_file_location,row.names=1)
results2<-read.csv("~/Desktop/Compare/LUAD_Neigh_45_3_TPM.matrix.csv_results_final.csv",row.names=1)
mat1<-as.matrix(fread(pii_values_file_location,data.table=FALSE))
mat2<-as.matrix(fread("~/Desktop/Compare/LUAD_Neigh_45_3_TPM.matrix.csv-181473-2015-06-19_pii_values_udi.csv",data.table=FALSE))

dim(mat2)
dim(mat1)  
genes_of_interest_1<<-columns_cutoff(results1,q_value_cutoff=.05,min_frac=-1,max_frac=1,equal=FALSE)
genes_of_interest_2<<-columns_cutoff(results2,q_value_cutoff=.05,min_frac=-1,max_frac=1,equal=FALSE)
print(paste("Numer of genes selected 1:",length(genes_of_interest_1)))
print(paste("Numer of genes selected 2:",length(genes_of_interest_2)))


#intersect(colnames(mat2),genes_of_interest_2)
#mat2[,intersect(colnames(mat2),genes_of_interest_2)]


jsd_mat_slow<-JSD_matrix(mat1,mat2,cores = 2,d = 20,fast =FALSE)
jsd_mat_fast<-JSD_matrix(mat1,mat2,cores = 2,d = 20,fast =TRUE)

#colnames(jsd_mat_fast)<-genes_of_interest_1[1:ncol(jsd_mat_fast)]
#rownames(jsd_mat_fast)<-genes_of_interest_2[1:nrow(jsd_mat_fast)]

print("Are slow and fast matrices the same?:")
#print(identical(as.numeric(jsd_mat_slow),as.numeric(jsd_mat_fast)))


JSD_results_file<-paste0("JSD_matrix-",Sys.Date(),".csv")
write.csv(jsd_mat_fast,JSD_results_file,row.names=TRUE) #Initializing rolling log file
#write.csv(jsd_mat_slow,JSD_results_file,row.names=TRUE) #Initializing rolling log file

colnames(jsd_mat_fast)<-NULL
rownames(jsd_mat_fast)<-NULL

colnames(jsd_mat_slow)<-NULL
rownames(jsd_mat_slow)<-NULL



