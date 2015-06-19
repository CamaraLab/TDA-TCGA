#setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/SC-TDA/UDI")
#Preping environment, loading necessary libraries
#debug(permute)

library(igraph)
library(rgexf)
library(jsonlite)
library(parallel)
library(getopt)
library(data.table)


#Setting defaults for debug mode
arg<-list("LUAD_Neigh_45_3","TPM.matrix.csv","1:10",200,detectCores(),TRUE,TRUE,4)
names(arg)<-c("name","matrix","columns","permutations","cores","log2","fdr","chunk")

#Argument section handling
spec = matrix(c(
  "name", "n", 1, "character",
  "matrix", "m",1,"character",
  "columns", "c",1,"character",
  "genes", "g",1,"character",
  "permutations","p",2,"integer",
  "cores","q",1,"integer",
  "log2","l",2,"logical",
  "fdr","f",2,"logical",
  "chunk","k",2,"integer"  
), byrow=TRUE, ncol=4)

#arg<-getopt(spec) #Conmment this line for debug mode

if ( is.null(arg$permutations ) ) {arg$permutations= 500}
if ( is.null(arg$log2 ) ) {arg$log2= TRUE}
if ( is.null(arg$fdr ) ) {arg$fdr= TRUE}
if ( is.null(arg$cores ) ) {arg$cores= detectCores()}
if ( is.null(arg$chunk ) ) {arg$chunk= 200}
if ( is.null(arg$columns ) ) {arg$columns= "all"}



#Loading matrix file to memory and log transforming if log2=TRUE
matrix_full_name<-arg$matrix
print (paste0("Loading ",matrix_full_name, " file to memory"))
matrix1<-fread(matrix_full_name,data.table=FALSE)
rownames(matrix1)<-matrix1[,1]
matrix1<-matrix1[,-1]
matrix1<-as.matrix(matrix1) #Converting to matrix from data.frame -> increases speed dramatically
if (arg$log2==TRUE) {matrix1<-(2^matrix1)-1} else matrix1<-as.matrix(matrix1)




column_range<-function(col_range)
  #Gets column range from arg$column and parse it.
{
  
  if (col_range=="all") #Check if all is supplied
  {
    x<-1:ncol(matrix1)
    
  } else if (grepl(":",col_range)==TRUE) { #If range x:x is supplied
    
    x<-as.numeric(strsplit(col_range,":")[[1]][1]:strsplit(col_range,":")[[1]][2])
    
  } else if (suppressWarnings(!is.na(as.numeric(col_range)))) { # If column is all numbers  
    x<-as.numeric(col_range)  
  } else if (!is.na(match(col_range,colnames(matrix1)))) { #If Gene name exist, convert to column number
    x<-which(colnames(matrix1)==col_range)       
  } else stop ("Column name not found")
  
  
  
  #Validating column in range:
  if (max(x)>ncol(matrix1)) 
  {
    stop (paste0("Columns out of range"," Available range: 1:",ncol(matrix1)))
  }
  return(x)
}

#Extracting columns from arguments
columns<-column_range(arg$columns)

#Parsing and loading, gexf(edge file) and json (nodes file) to memory.
print ("Parsing json and gexf files")
graph_name<-arg$name
json_file<-paste0(graph_name,".json")
gexf_file<-paste0(graph_name,".gexf")
graph_gexf<-read.gexf(gexf_file)
graph_igraph<-gexf.to.igraph(graph_gexf)

nodes<-fromJSON(json_file) #List of samples within nodes
names(nodes)<-as.numeric(names(nodes))+1 # Starting node is 1 - rownames
edges<-get.edgelist(graph_igraph,names=FALSE) # List of nodes and edges

#Subsetting for the largest connected subgraph
cluster_list<-clusters(graph_igraph)
largest_cluster_id<-which.max(cluster_list$csize)
largest_cluster_nodes<-which(cluster_list$membership==largest_cluster_id)
edges<-subset(edges,edges[,1] %in% (largest_cluster_nodes))
nodes<-nodes[largest_cluster_nodes]


#relabling nodes and updating edges accordingly
nodes_relabling_table<-cbind(largest_cluster_nodes,1:length(largest_cluster_nodes))
rownames(nodes_relabling_table)<-largest_cluster_nodes
colnames(nodes_relabling_table)<-c("original","new")
edges<-apply(edges,2,function(x) sapply(x,function(old_label) old_label<-nodes_relabling_table[as.character(old_label),"new"]))
names(nodes)<-sapply(names(nodes),function(old_label) old_label<-nodes_relabling_table[as.character(old_label),"new"])

#Relabling samples
nodes<-lapply(nodes,function(x) x+1) #Adding one to each smaple label min(sample)==1 
samples<-unique(unlist(nodes))
samples_relabling_table<-cbind(samples,1:length(samples)) #Creating relabling table
rownames(samples_relabling_table)<-samples
colnames(samples_relabling_table)<-c("original","new")
nodes<-lapply(nodes,function (x) sapply(x, function (old_sample) old_sample<-samples_relabling_table[as.character(old_sample),"new"])) #Updating samples according to relabling table

#Initializing rolling results file
p_table<-NULL
col_rolling<-t(c("Gene","c_score","p_value","pi_mean","pi_sd","pi_fraction"))
unique_id<-round(runif(1, min = 111111, max = 222222),0)
file_prefix<-paste0(arg$name,"_",arg$matrix,"-",unique_id,"-",Sys.Date())

print(paste("File unique identifier:",unique_id))
write.table(col_rolling,paste0(file_prefix,"_results_rolling.csv"),sep=",",col.names=FALSE,row.names=FALSE)




connectivity_pvalues_table<-function (column,permutations,cores)
  # Generates connectivity p_values tables. 
{
  
  
  
  #Spliting job into chunks of chunk_size columns
  split.column<-split(column,ceiling(seq_along(column)/chunk_size))
  
  sapply(split.column,function (chunk) 
    #calculating c-scores and p-values for each chunk of columns
  {
    print (paste0(Sys.time()," Calculating c-scores and p-values for rows: ",chunk[1],"-",chunk[chunk_size]))
    
    partial_table<-parSapply(cl,chunk,function(x) 
      #partial_table<-sapply(chunk,function(x)       
    {
      #Check for identica
      p_connectivity(nodes,x,permutations)
    })
    
    partial_table<-t(partial_table)
    rownames(partial_table)<-colnames(matrix1)[chunk]
    write.table(partial_table,paste0(file_prefix,"_results_rolling.csv"),append=TRUE,sep=",",col.names=FALSE)
    
  })
  
  complete_table<-read.csv(paste0(file_prefix,"_results_rolling.csv"),row.names=1,header=TRUE)
  return(complete_table)
  
}


p_connectivity<-function(nodes,column,permutations)
  #Given, nodes,column of interest, and n permutations, return c value and it's corresponding p-value
{
  if (sum(matrix1[,column]==0)==nrow(matrix1)) #If column is all zeros Skip all permutations 
  {p_table<-rep(0,length(col_rolling)-1)} else 
  {
    pii<-pii_calc(nodes,column) #Vector of Pi values for particular columns in every node
    pii_mean<-mean(pii)
    pii_sd<-sd(pii)
    pii_frac<-sum(pii!=0)/length(pii)
    c_score<-c_calc(edges,pii) # Connectivity value for specific column across the graph
    p_table<-c(c_score,p_value(column,permutations,c_score),pii_mean,pii_sd,pii_frac)
  }   
  return(p_table)
} 

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


p_value<-function (column,permutations,c_score) {
  #Retunrs p-value of a specific c-score of a specific column  
  if (permutations==0) {return(0)} else {
    c_vec<-NULL
    
    all_samples<-unlist(nodes,use.names=FALSE)
    samples<-unique(all_samples)
    perm_dic<-matrix(NA,length(samples),permutations) #rows=  cols=
    perm_dic<-apply(perm_dic,2,function(x) x<-sample(samples))
    permuted_samples<-matrix(NA,length(all_samples),permutations) # rows= #cols =
    
    
    
    for (i in seq_along(samples)) {
      samples_to_replace<-which (all_samples==i)
      for (j in seq_along(samples_to_replace))
        permuted_samples[samples_to_replace[j],]<-perm_dic[i,]   
    }
    
    for (i in 1:permutations)
    {
      permutation<-relist(permuted_samples[,i],skeleton=nodes)
      #permutation<-permute(nodes) #Permuted nodes 0.02
      pii_perm<-pii_calc(permutation,column) #Pi vector of specific column in permuted nodes 0.12
      c_perm<-c_calc(edges,pii_perm) #0.2
      c_vec[i]<-round(c_perm,5)
    }
    p<-length(which (c_vec>c_score))/length(c_vec)
    
  }
}






c_calc<-function(edges,pii)
  #Calculate connectivity value of a prticular column, 
  #pi values,based on current edges structure. (sim Adjacency matrix)
  #c=pi*pj*Aij()
{
  
  #microbenchmark(c1<-apply(edges,1,function (x) return (pii[x[1]]*pii[x[2]])))
  #c1<-apply(edges,1,function (x) return (pii[x[1]]*pii[x[2]]))
  c<-pii[edges[,1]]*pii[edges[,2]]
  
  #microbenchmark(c2<-pii[edges[,1]]*pii[edges[,2]])
  c<-sum(c)*2 #Multiple by 2 for 2 way edges.
  c<-c*length(largest_cluster_nodes)/(length(largest_cluster_nodes)-1) #Normalizing for graph size
  return(c)
  
}


permute<-function(nodes)
  #Permuting samples within nodes, maintaining skeleton and samples labels
{
  a<-unlist(nodes,use.names=FALSE)
  permuted<-sample(unique(a))
  d<-sapply(a,function (original) permuted[original])
  l<-relist(d,skeleton=nodes)
  return(l)
}



print(paste("Number of columns:",length(columns)))
print(paste("Number of permutations:",arg$permutations))


#parallel computing prep
print (paste("Acquiring",arg$cores,"cpu cores"))
cl <- makeCluster(as.numeric(arg$cores))
varlist=c("p_connectivity","arg","p_value","permute","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns","samples_relabling_table")
clusterExport(cl=cl, varlist=varlist,envir=environment())


ptm<<-proc.time() #Sytem time stamp for running calculation
chunk_size<-arg$chunk    #Size of column chunk


##########################################################
##########################################################
##########################################################
#########################################################

#This is the RUN command:

ans<-connectivity_pvalues_table(column=columns,permutations=arg$permutations,cores=arg$cores)
#ans<-as.matrix(ans) # Do not transform ans to matrix - in that case n=1 columns file is bad
non_zero_columns<-which(ans[,"c_score"]!=0) #Filtering for non-zero c-value columns
ans<-ans[non_zero_columns,] #removing nodes with zero c-score in final file


##########################################################
###########################################################
##########################################################
###########################################################
##########################################################
###########################################################


pii_values_table<-function(nodes,non_zero_columns) {
  #Returns pii_values of requested columns
  pii_values <-parSapply(cl,non_zero_columns,function (x) pii_calc (nodes,x))
  colnames(pii_values)<-colnames(matrix1)[non_zero_columns]
  rownames(pii_values)<-(1:length(nodes))
  return(pii_values)
}



print ("Correcting for multiple testings:")
#Writing final file if TRUE sorting by q-values, otherwise just sorting p-values
if (arg$fdr==TRUE)
{
  p_values<-ans[,2]
  q_values<-p.adjust(p_values,"fdr")
  ans<-cbind(ans,q_values)
  colnames(ans)[ncol(ans)]<-("q_value")
  ans<-ans[order(ans[,"q_value"]),]
} else  ans<-ans[order(ans[,"p_value"]),]


print ("Writing results_final.csv file to disk:")
write.csv(ans,paste0(file_prefix,"_results_final.csv"))



run_t<-proc.time()-ptm #Calculationg run time
print(paste("Runtime in seconds:",run_t[3]))
print(paste("Average p_value calc per column:",run_t[3]/length(columns)))
print(paste("Speed index (calc time for 500 permutations):",run_t[3]*500/arg$permutations/length(columns)))


print("Writing pii_values file:")
columns_of_interest<-match(rownames(ans),colnames(matrix1))
pii_values<-pii_values_table(nodes,columns_of_interest)
write.table(cbind(1:length(nodes),pii_values),paste0(file_prefix,"_pii_values.csv"),sep=",",row.names=FALSE)

print("Releasing cores")
stopCluster(cl) # Releasing acquired CPU cores
print(paste("File unique identifier:",unique_id))