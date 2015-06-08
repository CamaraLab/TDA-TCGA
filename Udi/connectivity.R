setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/SC-TDA/UDI")
#Preping environment, loading necessary libraries

library(igraph)
library(rgexf)
library(jsonlite)
library(parallel)
library(getopt)
  

#Setting defaults for debug mode
arg<-list("LUAD_Neigh_45_3","TPM.matrix.light.csv","10",20,detectCores(),TRUE,TRUE,200)
names(arg)<-c("name","matrix","columns","permutations","cores","log2","fdr")

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

arg<-getopt(spec)

if ( is.null(arg$permutations ) ) {arg$permutations= 500}
if ( is.null(arg$columns ) ) {arg$columns= 1}
if ( is.null(arg$log2 ) ) {arg$log2= TRUE}
if ( is.null(arg$fdr ) ) {arg$fdr= TRUE}
if ( is.null(arg$cores ) ) {arg$cores= detectCores()}
if ( is.null(arg$chunk ) ) {arg$chunk= 200}
if ( is.null(arg$columns ) ) {arg$columns= "all"}



#Loading matrix file to memory and log transforming if log2=TRUE
matrix_full_name<-arg$matrix
print (paste0("Loading ",matrix_full_name, " file to memory"))
matrix1<-read.csv(matrix_full_name,row.names=1)
if (arg$log2==TRUE) {matrix1<-(2^matrix1)-1}



#Extracting columns to analyze from arguments

#arg$columns<-"all"
#column<-column_range(arg$columns)
#column   
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

columns<-column_range(arg$columns)
#column   
    
##if (grepl(":",arg$columns)==TRUE) #If range x:x is supplied
  ##{  
  #columns<-as.numeric(strsplit(arg$columns,":")[[1]][1]:strsplit(arg$columns,":")[[1]][2])
#} else 
#    columns<-as.numeric(arg$columns)





#if gene name to column number
#if (!is.null(arg$genes)) 
#{
#  if (!is.na(match(arg$genes,colnames(matrix1)))) #Check if gene name exist
##    {columns<-which(colnames(matrix1)==arg$genes)} else stop("Column name not found")
#}

#}

#Parsing and loading, gexf(edge file) and json (nodes file) to memory.
print ("Parsing json and gexf files")
graph_name<-arg$name
json_file<-paste0(graph_name,".json")
gexf_file<-paste0(graph_name,".gexf")
graph_gexf<-read.gexf(gexf_file)
graph_igraph<-gexf.to.igraph(graph_gexf)

nodes<-fromJSON(json_file) #List of samples within nodes
edges<-get.edgelist(graph_igraph,names=FALSE) # List of nodes and edges

#Subsetting for the largest connected subgraph
cluster_list<-clusters(graph_igraph)
largest_cluster_id<-which.max(cluster_list$csize)
largest_cluster_nodes<-which(cluster_list$membership==largest_cluster_id)
edges<-subset(edges,edges[,1] %in% largest_cluster_nodes)



#Initializing rolling results file
p_table<-NULL
col1<-t(c("Gene","c-score","p-value"))
write.table(col1,paste0(arg$name,"_",arg$matrix,"_results_rolling.csv"),sep=",",col.names=FALSE,row.names=FALSE)




connectivity_pvalues_table<-function (column,permutations,cores)
  # Generates connectivity p_values tables. 
{
  
  #parallel computing prep
  print ("Acquiring cpu cores")
  cl <- makeCluster(as.numeric(cores))
  varlist=c("p_connectivity","arg","p_value","permute","c_vector","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes")
  clusterExport(cl=cl, varlist=varlist,envir=environment())
  
  ptm<<-proc.time() #Sytem time stamp for running calculation
  chunk_size<-arg$chunk    #Size of column chunk
  
  #Spliting job into chunks of chunk_size columns
  split.column<-split(column,ceiling(seq_along(column)/chunk_size))
  
  sapply(split.column,function (chunk) 
  #calculating c-scores and p-values for each chunk of columns
  {
    print (paste0(Sys.time()," Calculating C-scores and P-values for rows: ",chunk[1],"-",chunk[chunk_size]))
    partial_table<-parSapply(cl,chunk,function(x) 
          {
             p_connectivity(nodes,x,permutations)
          })
    
    partial_table<-t(partial_table)
    rownames(partial_table)<-colnames(matrix1)[chunk]
    p_table<<-rbind(p_table,partial_table)
    write.table(partial_table,paste0(arg$name,"_",arg$matrix,"_results_rolling.csv"),append=TRUE,sep=",",col.names=FALSE)
    
  })
  
  print("Releasing cores")
  
  stopCluster(cl) # Releasing acquired CPU cores
  return (p_table)
}


p_connectivity<-function(nodes,column,permutations)
  #Given, nodes,column of interest, and n permutations, return c value and it's p-value
{
  if (sum(matrix1[,column]==0)==nrow(matrix1)) #If column is all zeros Skip all permutations 
  {p_table<-c(0,0)} else 
  {
  pii<-pii_calc(nodes,column) #Pi values of particular columns
  c1<-c_calc(edges,pii) #Connectivity value (c) before permuting
  c1<-c1*length(largest_cluster_nodes)/(length(largest_cluster_nodes)-1) #Fine tunning c-value
  c<-c_vector(nodes,column,permutations) #Connectivity vector - all c values across all permutations
  p_table<-c(c1,p_value(c1,c))
  }   
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
  if (arg$log2==TRUE) {
    ei<-sapply(nodes,function (x) log2(1+mean(matrix1[x+1,column]))) 
        
  } else ei<-sapply(nodes,function (x) mean(matrix1[x+1,column]))
  
  total_ei<-sum(ei)  
  if (total_ei==0) {pii<-ei} else pii<-ei/total_ei
  #pii<-ei/total_ei
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
  #c=pi*pj*Aij()
{
  
  c<-apply(edges,1,function (x) pii[x[1]]*pii[x[2]])
  c<-sum(c)*2 #Multiple by 2 for 2 way edges.
  return(c)
  
}



permute<-function(nodes)
  #Permuting rows within nodes, maintaining skeleton and samples labels
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



##########################################################
##########################################################
##########################################################
#########################################################

#This is the RUN command:

ans<-connectivity_pvalues_table(column=columns,permutations=arg$permutations,cores=arg$cores)


##########################################################
###########################################################
##########################################################
###########################################################
##########################################################
###########################################################





print ("Correcting for multiple testings:")
#Writing final file if TRUE, otherwise just sorting p-values
if (arg$fdr==TRUE)
{
  p_values<-ans[,2]
  q_values<-p.adjust(p_values,"fdr")
  ans<-cbind(ans,q_values)
  colnames(ans)<-c("c-score","p-value","q-value")
  ans<-ans[order(ans[,"q-value"]),]
  } else 
      {
        colnames(ans)<-c("c-score","p-value")
        ans<-ans[order(ans[,"p-value"]),]
        
      }
  print ("Writing results_final.csv file to disk:")
  write.csv(ans,paste0(arg$name,"_",arg$matrix,"_results_final.csv"))



run_t<-proc.time()-ptm #Calculationg run time
print(paste("Runtime in seconds:",run_t[3]))
print(paste("Number of columns:",length(columns)))
print(paste("Number of permutations:",arg$permutations))
print(paste("Average p_value calc per column:",run_t[3]/length(columns)))
