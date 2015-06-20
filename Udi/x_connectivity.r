library(microbenchmark)
permutations<-500
columns<-5
nodes
all_samples<-unlist(nodes,use.names=FALSE)
samples<-unique(all_samples)



perm_dict<-as.list(rep(NA,columns))
dict<-matrix(NA,length(samples),permutations) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
perm_dict<-lapply(perm_dict,function(x) x<-apply(dict,2,function(x) x<-sample(samples)))


permuted_values<-matrix(NA,length(all_samples),permutations) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_dic
permuted_samples<-matrix(NA,length(all_samples),permutations) # rows= samples across graph,cols = permutation ID,flash= corresponding permuted sample ID from perm_di


translated_samples<-as.list(rep(NA,columns))
translated_values<-as.list(rep(NA,columns))
e_list<-as.list(rep(NA,columns))
pi_list<-as.list(rep(NA,columns))

translated_samples<-lapply(perm_dict,function(x) translate_samples(x,all_samples))
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

e_list<-lapply(translated_values,function(x) e_matrix(nodes,as.matrix(x)))


permuted_values<-as.matrix(translated_values[[1]])
e_matrix<-function(nodes,permuted_values) {
  e_matrix<-sapply(nodes,function (x) 
    apply(permuted_values[x,],2,mean))
  e_matrix<-t(e_matrix)
}
  
  
pi_list<-lapply(e_list,function(x) pii_matrix(as.matrix(x)))


pii_matrix<-function(e_matrix){
  pii_matrix<-apply(e_matrix,2,function (x)
    if (sum(x)==0) {pii_column<-x} else
      pii_column<-x/sum(x)
  )
}
  
    
    if (sum(x)==0) {pii<-ei} else pii<-ei/total_ei


pii_calc<-function(nodes,column)
  {  #Calculate pi value(mean) for a particular column in a particular node.
    # Input: nodes list and one specific column(i) of interest, output: pi
    if (arg$log2==TRUE) {
      ei<-sapply(nodes,function (x) log2(1+mean(matrix1[samples_relabling_table[x,1],column]))) 
      
    } else ei<-sapply(nodes,function (x) mean(matrix1[samples_relabling_table[x,1],column]))
    
    total_ei<-sum(ei)
    if (total_ei==0) {pii<-ei} else pii<-ei/total_ei
    return(pii)
  }


