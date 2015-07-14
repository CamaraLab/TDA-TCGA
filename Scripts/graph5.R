#setwd(".")
#setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
library(igraph)
library(rgexf)
library(jsonlite)
library(parallel)

connectivity_pvalues_table<-function (graph_name,matrix_full_name,column,permutations,FDR,log_scale,cores)
  # Generates connectivity p_values tables. 
{
  n_permutations<<-permutations #Used for console printing only
  n_columns<<-length(column) #Used for console printing only
  print (paste0("[1] Loading ",matrix_full_name," file to memory"))
  matrix1<<-read.csv(matrix_full_name,row.names=1)
  #matrix1<<-matrix1[,1:500]
  if (log_scale) {matrix1<<-(2^matrix1)-1} 
  
  print ("[2] Parsing json and gexf files")
  json_file<-paste0(graph_name,".json")
  gexf_file<-paste0(graph_name,".gexf")
  graph_gexf<-read.gexf(gexf_file)
  graph_igraph<-gexf.to.igraph(graph_gexf)
  
  nodes<-fromJSON(json_file) #List of samples within nodes
  edges<-get.edgelist(graph_igraph,names=FALSE) # List of nodes and edges
  
  #parallel computing prep
  print ("Acquiring cpu cores")
  cl <- makeCluster(cores)
  varlist=c("permutations","column","parSapply","p_value","permute","c_vector","edges","nodes","matrix1","p_connectivity","cl","pii_calc","c_calc")
  clusterExport(cl=cl, varlist=varlist,envir=environment())
  
  ptm<<-proc.time()
  chunk_size<-50
  #Spliting job into chunks of 20(chunk) columns
  split.column<-split(column,ceiling(seq_along(column)/chunk_size))
  
  sapply(split.column,function (chunk) 
  #p_table<-parSapply(cl,split.column,function (chunk) 
  {
    print (paste0(Sys.time()," Calculating C-scores and P-values for rows: ",chunk[1],"-",chunk[chunk_size]))
    partial_table<-parSapply(cl,chunk,function(x) 
    #partial_table<-sapply(chunk,function(x)   
          {
             p_connectivity(nodes,x,permutations)
          })
    
    #print (paste0(Sys.time(),"Updating results.csv for rows: ",chunk[1],"-",chunk[20]))
    partial_table<-t(partial_table)
    rownames(partial_table)<-colnames(matrix1)[chunk]
    p_table<<-rbind(p_table,partial_table)
    write.table(partial_table,"results.csv",append=TRUE,sep=",",col.names=FALSE)
    
  })
  
  
  print("Releasing cores")
  stopCluster(cl)
  
  return (p_table)
}


p_connectivity<-function(nodes,column,permutations)
  #Given, nodes,column of interest, and n permutations, return c value and it's p-value
{
  pii<-pii_calc(nodes,column) #Pi values of particular columns
  c1<-c_calc(edges,pii) #Connectivity value (c) before permuting
  c<-c_vector(nodes,column,permutations) #Connectivity vector - all c values across all permutations
  p_table<-c(c1,p_value(c1,c))
  return(p_table)
} 

c_vector<-function(nodes,column,permutations)
  #Generate vector of connectivity values each element is connectivity value for each permutation
{
  c_vec<-NULL
  for (i in 1:permutations)
  {
    permutation<-permute(nodes) #Permuted nodes 0.02
    pii_perm<-pii_calc(permutation,column) #Pi vector of specific column in permuted nodes 0.12
    c<-c_calc(edges,pii_perm) #0.2
    c_vec[i]<-round(c,5)
  }
  return(c_vec)
  
}

pii_calc<-function(nodes,column)
{  #Calculate pi value(mean) for a particular column in a particular node.
  # Input: nodes list and one specific column(i) of interest, output: pi
  #if (log_scale==TRUE) 
    ei<-sapply(nodes,function (x) log2(1+mean(matrix1[x+1,column])))
  #else
   # ei<-sapply(nodes,function (x) mean(TPM.matrix[x+1,column]))
  #ei<-parSapply(cl,nodes,function (x) mean(TPM.matrix[x+1,column]))
  total_ei<-sum(ei)  
  pii<-ei/total_ei
  return(pii)
}



p_value<-function (c1,c_vec)
{
  # Given pi value and connectivity vector it returns p-value for the pi
  p<-length(which (c_vec>c1))/length(c_vec)
}

c_calc<-function(edges,pii)
  #Calculate connectivity value of a prticular column, 
  #pi values,based on current edges structure. (sim Adjacency matrix)
  #c=pi*pj*Aij
{
  c<-apply(edges,1,function (x) pii[x[1]]*pii[x[2]])
  con<-sum(c)*2 #Multiple by 2 for 2 way edges.
}



permute<-function(nodes)
  #Permuting rows within nodes, maintaining skeleton
{
  a<-unlist(nodes,use.names=FALSE)
  original<-unique(a)
  permuted<-sample(original)
  conversion<-cbind(original,permuted)
  conversion<-conversion[order(conversion[,"original"]),]
  
  d<-sapply(a,function (original)
  {
    a.value<-conversion[original+1,"permuted"]
  })
  
  l<-relist(d,skeleton=nodes)
  return(l)
}



p_table<<-NULL
col1<<-t(c("Row","c-score","p-value","q-value"))
write.table(col1,"results.csv",sep=",",col.names=FALSE,row.names=FALSE)

##########################################################
##########################################################
##########################################################
#########################################################

#This is the RUN command:

ans<-connectivity_pvalues_table("LUAD_Neigh_45_3","TPM.matrix.light.csv",column=1:5,permutations=10,FDR=TRUE,log_scale=TRUE,cores=4)
#ans<-connectivity_pvalues_table("LUAD_MDS_30_3","TPM.matrix.csv",column=1:2,permutations=100,FDR=TRUE,cores=4)

##########################################################
###########################################################
##########################################################
###########################################################
##########################################################
###########################################################


run_t<-proc.time()-ptm #Calculationg run time


print ("Correcting for multiple testings:")

p_values<-ans[,2]
q_values<-p.adjust(p_values,"fdr")
ans<-cbind(ans,q_values)
colnames(ans)<-c("c-score","p-value","q-value")
ans<-ans[order(ans[,"q-value"]),]


print ("[4] Writing results_FDR.csv file to disk")
write.csv(ans,"results_FDR.csv")

print(paste("Runtime in seconds:",run_t[3]))
print(paste("Number of columns:",n_columns))
print(paste("Number of permutations:",n_permutations))
print(paste("Average p_value calc per column:",run_t[3]/n_columns))
