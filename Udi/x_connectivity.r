library(microbenchmark)
permutations<-500
columns<-20
nodes
all_samples<-unlist(nodes,use.names=FALSE)
samples<-unique(all_samples)
edges1<-edges[,1]
edges2<-edges[,2]
num_nodes<-length(nodes)
matrix2<-matrix1[,7]

dict<-matrix(NA,length(samples),permutations+1) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
#dict[,1]<-seq_along(samples) #
perm_dict<-as.list(rep(NA,columns))
perm_dict<-lapply(perm_dict,function(x) {
    x<-apply(dict,2,function(x) x<-sample(samples))
    x[,1]<-1:length(samples)
    return(x)
  }) 



permuted_values<-matrix(NA,length(all_samples),permutations+1) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_dic
#permuted_samples<-matrix(NA,length(all_samples),permutations) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_di


#translated_samples<-as.list(rep(NA,columns))
#translated_values<-as.list(rep(NA,columns))
##e_list<-as.list(rep(NA,columns))
#pi_list<-as.list(rep(NA,columns))
#c_list<-as.list(rep(NA,columns))

#translated_samples<-lapply(perm_dict,function(x) translate_samples(x,all_samples))
column<-0
system.time(translated_values<-lapply(perm_dict,function(x) {
  column<-column+1
  translate_values(x,all_samples,column)
  #translate_values(x,all_samples)

}))

#microbenchmark(matrix2<-matrix1[,1])
#x2<-translate_samples(perm_dict[[2]],all_samples)


 


cl <- makeCluster(as.numeric(arg$cores))
varlist=c("e_matrix1","p_connectivity","arg","p_value","permute","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns","samples_relabling_table")
clusterExport(cl=cl, varlist=varlist,envir=environment())


e_list1<-lapply(translated_values,function(x) e_matrix1(nodes,as.matrix(x)))
e_list1<-parLapply(cl,translated_values,function(x) e_matrix1(nodes,as.matrix(x)))
#system.time(e_list1<-parLapply(cl,translated_values,function(x) e_matrix1(nodes,as.matrix(x))))




#colMeans(permuted_values,x)))
  
pi_list<-lapply(e_list1,function(x) pii_matrix(as.matrix(x)))
system.time(c_vec_list<-lapply(pi_list,function (pi_matrix) c_calc(as.matrix(pi_matrix)))) # return a list where each element contains C_Vector for a column


# Go over each pii_matrix columns 

c_scores<-sapply(c_vec_list,function (c_vec) c_vec[1])

length(c_scores)

#hugh_translated<-matrix(NA,length(all_samples),columns*permutations)
#hugh_translated<-NULL
#system.time(for (i in seq_along(translated_values)) {
#  hugh_translated<-cbind(hugh_translated,translated_values[[i]])
#})