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
arg<-list("LUAD_Neigh_45_3","TPM.matrix.light.csv","1",10,detectCores(),FALSE,TRUE,4)
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



perm_values<-function(dict_matrix,samples_relabling_table,column) {
  #Takes dictionary matrix and all samples- returns translated_matrix with corersponding values
  #z<-apply(y,2,function (x) x<-matrix2[samples_relabling_table[x,1]])
  y<-apply(dict_matrix,2,function(perm_column) perm_column<-matrix1[samples_relabling_table[perm_column,1],column])
  #y<-apply(dict_matrix,2,function(perm_column) perm_column<-matrix2[samples_relabling_table[perm_column,1]])
  
}
  

e_matrix<-function(nodes,translated_values) {
  e_matrix<-sapply(nodes,function (x) {
    s<-0
    for (i in 1:length(x)) {
      #s<-s+translated_values[x,][i,]

      s<-s+translated_values[x,][i,]
    }
    #e<-log2(1+s/length(x))
    e<-s/length(x)
  })
  e_matrix<-t(e_matrix) #Transposing to keep permutations as columns
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
permutations<-10
columns<-10
#nodes
all_samples<-unlist(nodes,use.names=FALSE)
samples<-unique(all_samples)
edges1<-edges[,1]
edges2<-edges[,2]
num_nodes<-length(nodes)
matrix2<-matrix1[,8]
column<-8

dict<-matrix(NA,length(samples),permutations+1) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
#dict[,1]<-seq_along(samples) #
perm_dict_list<-as.list(rep(NA,columns)) #Creating list of "column" elements
perm_dict_list<-lapply(perm_dict_list,function(x) { #Each element in perm_dict list will get permutation matrix 
  x<-apply(dict,2,function(x) x<-sample(samples))
  x[,1]<-1:length(samples)
  return(x)
}) 


#system.time(perm_values_list<-lapply(perm_dict_list,function (x) x<-perm_values(x,samples_relabling_table))) #Translate permuted sample to corresponding value in matrix2


#perm_values_list<-lapply(perm_dict_list,function (x) x<-perm_values(x,samples_relabling_table))

perm_values_list<-as.list(rep(NA,columns)) #Creating list of "column" elements

for (column in 1:columns) {
  perm_values_list[[column]]<-perm_values(perm_dict_list[[column]],samples_relabling_table,column)
}

#perm_values_list<-lapply(perm_dict_list,function (x) x<-perm_values(x,samples_relabling_table,column))


#cl <- makeCluster(as.numeric(arg$cores))
#varlist=c("e_matrix","p_connectivity","arg","p_value","permute","edges","nodes","matrix1","pii_calc","c_calc","largest_cluster_nodes","col_rolling","columns","samples_relabling_table")
#clusterExport(cl=cl, varlist=varlist,envir=environment())


#system.time(e_list1<-lapply(perm_values_list,function(x) e_matrix(nodes,as.matrix(x))))
e_list1<-lapply(perm_values_list,function(x) e_matrix(nodes,as.matrix(x)))
#e_list1<-parLapply(cl,translated_values,function(x) e_matrix(nodes,as.matrix(x)))
#system.time(e_list1<-parLapply(cl,translated_values,function(x) e_matrix(nodes,as.matrix(x))))




#system.time(pi_list<-lapply(e_list1,function(x) pii_matrix(as.matrix(x))))
#system.time(c_vec_list<-lapply(pi_list,function (pi_matrix) c_calc_fast(as.matrix(pi_matrix)))) # return a list where each element contains C_Vector for a column
pi_list<-lapply(e_list1,function(x) pii_matrix(as.matrix(x)))
c_vec_list<-lapply(pi_list,function (pi_matrix) c_calc_fast(as.matrix(pi_matrix)))

# Go over each pii_matrix columns 

#system.time(c_scores<-sapply(c_vec_list,function (c_vec) c_vec[1]))
c_scores<-sapply(c_vec_list,function (c_vec) c_vec[1])

print(c_scores)
length(c_scores)




