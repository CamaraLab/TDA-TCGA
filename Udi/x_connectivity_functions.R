translate_values<-function(dict_matrix,all_samples) {
  #Takes dictionary matrix and all samples- returns translated_matrix with corersponding values
  for (i in samples) {
    samples_to_replace<-which (all_samples==i)
    for (j in samples_to_replace) {
      #permuted_samples[j,]<-perm_dict[i,]   
      permuted_values[samples_to_replace[j],]<-matrix2[dict_matrix[i,]]   
    
      #permuted_values[j,]<-matrix1[perm_dict[i,],column]
    }
  }
  return<-permuted_values
}


e_matrix1<-function(nodes,translated_values) {
  e_matrix1<-sapply(nodes,function (x) {
    s<-0
    for (i in 1:length(x)) {
      s<-s+translated_values[x,][i,]
    }
    #e<-log2(1+s/length(x))
    e<-s/length(x)
  })
  e_matrix1<-t(e_matrix1) #Transposing to keep permutations as columns
}
  



e_matrix<-function(nodes,translated_values) {
  e_matrix<-sapply(nodes,function (x) 
    apply(translated_values[x,],2,mean))
  e_matrix<-t(e_matrix)
}



pii_matrix<-function(e_matrix){ #Gets e_matrix retuens pii
  pii_matrix<-apply(e_matrix,2,function (e_column)
    if (sum(e_column)==0) {pii_column<-e_column} else
      pii_column<-e_column/sum(e_column)
  )
}




c_calc_fast<-function(pi_matrix)
  #Calculate connectivity value of a prticular column, 
  #pi values,based on current edges structure. (sim Adjacency matrix)
  #c=pi*pj*Aij()
{
  c_vector<-apply(pi_matrix,2,function(pi_column) 
                        c<-2*sum(pi_column[edges1]*pi_column[edges2]))
  
  c_vector<-c_vector*(num_nodes/(num_nodes-1))
} 
  
library(microbenchmark)
permutations<-20
columns<-1
#nodes
all_samples<-unlist(nodes,use.names=FALSE)
samples<-unique(all_samples)
edges1<-edges[,1]
edges2<-edges[,2]
num_nodes<-length(nodes)
matrix2<-matrix1[,8]

dict<-matrix(NA,length(samples),permutations+1) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
#dict[,1]<-seq_along(samples) #
perm_dict<-as.list(rep(NA,columns))
perm_dict<-lapply(perm_dict,function(x) {
  x<-apply(dict,2,function(x) x<-sample(samples))
  x[,1]<-1:length(samples)
  return(x)
}) 



permuted_values<-matrix(NA,length(all_samples),permutations+1) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_dic
permuted_samples<-matrix(NA,length(all_samples),permutations+1) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_di
permuted_samples[,1]<-1:length(all_samples)

for (i in seq_along(samples)) {
  samples_to_replace<-which (all_samples==i)
  for (j in samples_to_replace)
    permuted_samples[all_samples[j],]<-perm_dict[[1]][i,]   
}



#translated_samples<-as.list(rep(NA,columns))
#translated_values<-as.list(rep(NA,columns))
##e_list<-as.list(rep(NA,columns))
#pi_list<-as.list(rep(NA,columns))
#c_list<-as.list(rep(NA,columns))

#translated_samples<-lapply(perm_dict,function(x) translate_samples(x,all_samples))
#column<-0
translated_values<-lapply(perm_dict,function(x) {
  #column<-column+1
  #translate_values(x,all_samples,column)
  translate_values(as.matrix(x),all_samples)
    
})

#microbenchmark(matrix2<-matrix1[,1])
#x2<-translate_samples(perm_dict[[2]],all_samples)





#cl <- makeCluster(as.numeric(arg$cores))
#varlist=c("e_matrix1","p_connectivity","arg","p_value","permute","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns","samples_relabling_table")
#clusterExport(cl=cl, varlist=varlist,envir=environment())


e_list1<-lapply(translated_values,function(x) e_matrix(nodes,as.matrix(x)))
#e_list1<-parLapply(cl,translated_values,function(x) e_matrix1(nodes,as.matrix(x)))
#system.time(e_list1<-parLapply(cl,translated_values,function(x) e_matrix1(nodes,as.matrix(x))))




#colMeans(permuted_values,x)))

pi_list<-lapply(e_list1,function(x) pii_matrix(as.matrix(x)))
system.time(c_vec_list<-lapply(pi_list,function (pi_matrix) c_calc_fast(as.matrix(pi_matrix)))) # return a list where each element contains C_Vector for a column


# Go over each pii_matrix columns 

c_scores<-sapply(c_vec_list,function (c_vec) c_vec[1])
print(c_scores)
length(c_scores)

#hugh_translated<-matrix(NA,length(all_samples),columns*permutations)
#hugh_translated<-NULL
#system.time(for (i in seq_along(translated_values)) {
#  hugh_translated<-cbind(hugh_translated,translated_values[[i]])
#})


