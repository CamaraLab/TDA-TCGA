library(microbenchmark)
library (doParallel)
library(foreach)
permutations<-500
columns<-100
nodes
all_samples<-unlist(nodes,use.names=FALSE)
samples<-unique(all_samples)
edges1<-edges[,1]
edges2<-edges[,2]
num_nodes<-length(nodes)

perm_dict<-as.list(rep(NA,columns))
dict<-matrix(NA,length(samples),permutations) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
perm_dict<-lapply(perm_dict,function(x) x<-apply(dict,2,function(x) x<-sample(samples)))
matrix2<-matrix1[,7]

permuted_values<-matrix(NA,length(all_samples),permutations) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_dic
#permuted_samples<-matrix(NA,length(all_samples),permutations) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_di


#translated_samples<-as.list(rep(NA,columns))
translated_values<-as.list(rep(NA,columns))
e_list<-as.list(rep(NA,columns))
pi_list<-as.list(rep(NA,columns))
c_list<-as.list(rep(NA,columns))

#translated_samples<-lapply(perm_dict,function(x) translate_samples(x,all_samples))
translated_values<-lapply(perm_dict,function(x) translate_values(x,all_samples))
#x2<-translate_samples(perm_dict[[2]],all_samples)



translate_samples<-function(perm_dict,all_samples) {
#Takes dictionary matrix and all samples- returns translated_matrix with corersponding values
  for (i in samples) {
    samples_to_replace<-which (all_samples==i)
    for (j in samples_to_replace) {
      permuted_samples[j,]<-perm_dict[i,]   
      #permuted_values[samples_to_replace[j],]<-matrix2[perm_dic[i,]]   
    }
  }
  return<-permuted_samples
}
 
translate_values<-function(perm_dict,all_samples) {
  #Takes dictionary matrix and all samples- returns translated_matrix with corersponding values
  for (i in samples) {
    samples_to_replace<-which (all_samples==i)
    for (j in samples_to_replace) {
      #permuted_samples[j,]<-perm_dict[i,]   
      permuted_values[j,]<-matrix2[perm_dict[i,]]   
    }
  }
  return<-permuted_values
}

cl <- makeCluster(as.numeric(arg$cores))
varlist=c("e_matrix1","p_connectivity","arg","p_value","permute","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns","samples_relabling_table")
clusterExport(cl=cl, varlist=varlist,envir=environment())


e_list<-lapply(translated_values,function(x) e_matrix1(nodes,as.matrix(x)))
e_list1<-papply(cl,translated_values,function(x) e_matrix1(nodes,as.matrix(x)))
system.time(e_list1<-lapply(translated_values,function(x) e_matrix1(nodes,as.matrix(x))))

dim(translated_values[[1]])

length(translated_values)

#hugh_translated<-matrix(NA,length(all_samples),columns*permutations)
hugh_translated<-NULL
system.time(for (i in seq_along(translated_values)) {
  hugh_translated<-cbind(hugh_translated,translated_values[[i]])
})

#foreach (i=seq_along(translated_values)) %do%
#  hugh_translated<-cbind(hugh_translated,translated_values[[i]])
#})

m<-matrix(1:100,10,10)
m

e_matrix2<-function(nodes,hugh_translated) {
  e_matrix1<-sapply(c(2,3),function (x)
    apply(m[c(2,3),],1,sum)/length(2)

(m[c(2,3),]

e_matrix1<-function(nodes,hugh_translated)
  e_matrix1<-sapply(nodes,function (x) {
    s<-0
    for (i in 1:length(x)) {
      s<-s+hugh_translated[x,][i,]
    }
    e<-s/length(x)
  })

    e<-apply(hugh_translated[x,],1,function(y) {
      s<-s+y
      e<-s/length(x)
      return(e)
    })
    
  })

    

e_matrix<-function(nodes,hugh_translated) {
  e_matrix<-sapply(nodes,function (x) 
    apply(hugh_translated[x,],2,mean))
  #e_matrix<-t(e_matrix)
}
#colMeans(permuted_values,x)))
  
pi_list<-lapply(e_list,function(x) pii_matrix(as.matrix(x)))


pii_matrix<-function(e_matrix){ #Gets e_matrix retuens pii
  pii_matrix<-apply(e_matrix,2,function (x)
    if (sum(x)==0) {pii_column<-x} else
      pii_column<-x/sum(x)
  )
}
 

c_list<-lapply(pii_matrix,function (x) c_calc(as.matrix(x)))

# Go over each pii_matrix columns 

if (sum(x)==0) {pii<-ei} else pii<-ei/total_ei


c_calc<-function(pii)
  #Calculate connectivity value of a prticular column, 
  #pi values,based on current edges structure. (sim Adjacency matrix)
  #c=pi*pj*Aij()
{
  c<-pii[edges1]*pii[edges2]
  c<-sum(c)*2 #Multiple by 2 for 2 way edges.
  c<-c*num_nodes/(num_nodes-1) #Normalizing for graph size
  return(c)
  
}

