
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

JSD_matrix<-function (pii_values,){
  
  
}
  

ans<-parLapply(cl,seq_along(pool),function (x) {
  
  a1<-matrix(0,d,length(pool[[x]]))
  
  for (j in min(pool[[x]]):max(pool[[x]]))
    for (i in (1:d)) {
      if (i>=j) {
        cols<-j-(x-1)*length(pool[[1]])
        a1[i,cols]<-JSD2(pii1[,i],pii1[,j])
        #a1[cols,i]<-a1[i,cols]}
        
      }
    }
  #return(a1)
  return(a1)
})








P<-pii1[,1]
Q<-pii1[,2]




#Generating pii_file
cl <- makeCluster(as.numeric(arg$cores))
varlist=c("p_connectivity","arg","p_value","permute","c_vector","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns")
clusterExport(cl=cl, varlist=varlist,envir=environment())


pii1<-as.matrix(fread("/Users/uer2102/Desktop/Merging2/LUAD_Neigh_45_3_TPM.matrix.csv_pii_values.csv",data.table=FALSE))
pii1<-pii1[,-1]
pii2<-as.matrix(fread("/Users/uer2102/Desktop/Merging2/LUAD_Neigh_45_3_TPM.matrix.csv-192553-2015-06-17_pii_values.csv",data.table=FALSE))
pii2<-pii2[,-1]
identical(round(pii1,4),round(pii2,4))
identical(pii1,pii2)

pii_val<-as.matrix(fread("/Users/uer2102/Desktop/Merging2/LUAD_Neigh_45_3_TPM.matrix.csv_pii_values.csv",data.table=FALSE))
pii_val1<-as.matrix(fread("/Users/uer2102/Desktop/Merging2/LUAD_Neigh_45_3_TPM.matrix.csv-192553-2015-06-17_pii_values.csv",data.table=FALSE))
results<-read.csv("/Users/uer2102/Desktop/Merging2/LUAD_Neigh_45_3_TPM.matrix.csv_results_final-3.csv",row.names=1)
pii_val<-pii_val[,-1]
pii_val1<-pii_val1[,-1]
pii1<-pii[,1]

#View(results)
#sum(results$pi_fraction==1)
results<-results[results$pi_fraction!=1,]
results<-results[results$q_value<=.05,]
genes_of_interest<-rownames(results)
pii1<-pii1[,genes_of_interest]
dim(pii1)



columns_of_interest<-match(rownames(combined),colnames(pii_val))
pii<-pii_val[,columns_of_interest]

d<-10


a<-matrix(0,d,d)


pool<-d
pool<-clusterSplit(cl,1:d)
pool


#parSapply(cl,pool,function(j) {
sapply(pool,function(j) {
  
  for (i in 1:length(j)) 
        sapply (j,function (i) { a[i,j]<-JSD1(pii1[,i],pii1[,j]) })

  })
  

for (1 in min(x):(max(x)-1)) 
    for (j in (i+1):d) {
      a[i,j]<-JSD1(pii1[,i],pii1[,j])
      a[j,i]<-a[i,j]
      
  
})

}
    

d<-100
a<-matrix(0,d,d)
system.time(
for (i in 1:(d-1)) 
  for (j in (i+1):d) {
    a[i,j]<-JSD2(pii1[,i],pii1[,j])
    a[j,i]<-a[i,j]
  }
)

a_slow<-a
a_fast<-a

identical(a_slow,a_fast)
a_slow[,1]
a_fast[,1]
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

number_of_identical_columns<-sum(sapply(mylist,function(x) length(x)-1))

1000-number_of_identical_columns

sum(sapply(mylist,function(x) if (length(x)-1==0) {return(0)} else return(1)))



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
