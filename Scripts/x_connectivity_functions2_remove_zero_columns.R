#setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/SC-TDA/UDI")
#Preping environment, loading necessary libraries

library(igraph)
library(rgexf)
library(jsonlite)
library(parallel)
library(getopt)
library(data.table)


#Setting defaults for debug mode
arg<-list("LUAD_Neigh_45_3","TPM.matrix.light.csv","50",250,detectCores(),TRUE,TRUE,200)
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

arg<-getopt(spec) #Conmment this line for debug mode

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

#Extracting columns from arguments and removign zero columns from matrix
columns<-column_range(arg$columns)
zero_columns<-which(apply(matrix1,2,function (x) sum(x==0)==nrow(matrix1))) #If column is all zeros removing from later calcuations
columns<-setdiff(columns,zero_columns)


#matrix2<-matrix1[,non_zero_columns]


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
samples<-unique(unlist(nodes))
matrix1<-matrix1[samples_relabling_table[,1],] #Subseting matrix to contain only samples in first connected graph 


#Initializing rolling results file
unique_id<-round(runif(1, min = 111111, max = 222222),0)
file_prefix<-paste0(arg$name,"_",arg$matrix,"-",unique_id,"-",Sys.Date())
print(paste("File unique identifier:",unique_id))

#Info_cols is used to set inforation columns in output file as well as names for the variables that constitutes those columns
info_cols<-t(c("Genes","c_scores","p_values","pi_frac","n_samples","e_mean","e_sd")) 
write.table(info_cols,paste0(file_prefix,"_results_rolling.csv"),sep=",",col.names=FALSE,row.names=FALSE)

print("Writing out file with zero pi-fraction columns:")
write.table(names(zero_columns),paste0(file_prefix,"_results_zero_frac.csv"),sep=",",col.names=FALSE,row.names=FALSE)



#perm_values<-function(dict_matrix,samples_relabling_table,column,matrix) {
perm_values<-function(dict_matrix,column,matrix) {
  #Takes dictionary matrix and all samples- returns translated_matrix with corersponding values
  
  #y<-apply(dict_matrix,2,function(perm_column) perm_column<-matrix[samples_relabling_table[perm_column,1],column])
  y<-apply(dict_matrix,2,function(perm_column) perm_column<-matrix[perm_column,column])
  
}
  

e_matrix<-function(nodes,translated_values) {
  e_matrix<-sapply(nodes,function (x) {
    s<-0
    for (i in 1:length(x)) {
      s<-s+translated_values[x,][i,]
    }
    if (arg$log2==TRUE) {
      e<-log2(1+s/length(x))} else e<-s/length(x)
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
  
permutations<-arg$permutations
edges1<-edges[,1] #Nodes i
edges2<-edges[,2] #Nodes j
num_nodes<-length(nodes)

chunk_size<-arg$chunk #How many columns each CPU will be assigned at once
split.column<-split(columns,ceiling(seq_along(columns)/chunk_size))


print (paste("Preparing parallel environment. Acquiring",arg$cores,"Cores"))
cl <- makeCluster(as.numeric(arg$cores))
varlist=c("file_prefix","c_calc_fast","c_calc_fast","pii_matrix","e_matrix","edges1","edges2","samples","permutations","num_nodes","chunk_size","perm_values","arg","nodes","matrix1","largest_cluster_nodes","info_cols","columns","samples_relabling_table")
clusterExport(cl=cl, varlist=varlist,envir=environment())

print (Sys.time())
ptm<-proc.time()

ans<-parLapply(cl,split.column,function (columns_range)  {
  #calculating c-scores and p-values for each chunk of columns
    
  #REMOVE ZERO COLUMNS HERE matrix2<-matrix1[,columns]
  
  matrix2<-matrix1[,columns_range] #Subsetting matrix2 to have only non_Zero_columns as passed by columns_range
  dict<-matrix(NA,length(samples),permutations+1) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
  perm_dict_list<-as.list(rep(NA,length(columns_range))) #Creating list of "column" elements
  perm_values_list<-as.list(rep(NA,length(columns_range))) #Creating list of "column" elements
  
  perm_dict_list<-lapply(perm_dict_list,function(x) { #Each element in perm_dict list will get permutation matrix 
      x<-apply(dict,2,function(x) x<-sample(samples))
      x[,1]<-1:length(samples)
      return(x)
    }) 


    
    
    
    for (column in seq_along(columns_range)) {
      #perm_values_list[[column]]<-perm_values(perm_dict_list[[column]],samples_relabling_table,column,matrix2)
      perm_values_list[[column]]<-perm_values(perm_dict_list[[column]],column,matrix2)
    }
    
    
    e_list<-lapply(perm_values_list,function(x) e_matrix(nodes,as.matrix(x))) #columns_range elements in the list. Each element is a matrix representing pi_values of a gene. rows are nodes, columns are permutations. first column is non permuted.
    pi_list<-lapply(e_list,function(x) pii_matrix(as.matrix(x))) #columns_range elements in the list. Each element is a matrix representing pi_values of a gene. rows are nodes, columns are permutations. first column is non permuted.
    c_vec_list<-lapply(pi_list,function (pi_matrix) c_calc_fast(as.matrix(pi_matrix)))
    
    
    #e_values<-sapply(e_list,function (x) x[,1]))
    
    e_mean<-sapply(e_list,function (x) mean(x[,1])) #Taking mean of the first column (not permutations)
    e_sd<-sapply(e_list,function (x) sd(x[,1]))
    pi_values<-sapply(pi_list,function (x) x[,1]) #Is a matrix,each row is a node, each column in the matrix is pi values of a gene across nodes.
    pi_frac<-apply(pi_values,2,function (x) sum(x!=0)/length(x))
    n_samples<-apply(matrix2,2,function (x) sum(x!=0))
    c_scores<-sapply(c_vec_list,function (c_vec) c_vec[1])
    p_values<-sapply(c_vec_list,function(c_vec) {
    p_value<-sum(c_vec>c_vec[1])/permutations})
    Genes<-colnames(matrix2)  
  
  output<-cbind(Genes,c_scores,p_values,pi_frac,n_samples,e_mean,e_sd) #The variable names should match info_cols
  write.table(output,paste0(file_prefix,"_results_rolling.csv"),append=TRUE,sep=",",col.names=FALSE,row.names=FALSE)
  
  #Ans is columns_range length list. Each element contain to variables.
  #First variabl is "output" which is results matrix. the second element is pii_values matrix, its columns are genes and rows are nodes 
  ans<-list(output[],rbind(Genes,pi_values))
  return(ans) 
  
  })


#Generating results table. 
final_results<-NULL
for (i in 1:length(ans))
  final_results<-rbind(final_results,ans[[i]][[1]]) #Extracting output from ans


q_value<-p.adjust(final_results[,"p_values"],"fdr")
final_results<-cbind(final_results,q_value)
final_results<-final_results[order(final_results[,"q_value"]),]



#Generating pii_values table
pi_values_table<-NULL
for (i in 1:length(ans))
  pi_values_table<-cbind(pi_values_table,ans[[i]][[2]]) #Extracting pii_values from ans
colnames(pi_values_table)<-pi_values_table[1,]


#Writing final results and pii_values files
write.table(final_results,paste0(file_prefix,"_results_final.csv"),row.names=FALSE,sep=",")
write.table(pi_values_table,paste0(file_prefix,"_pii_values.csv"),sep=",",row.names=FALSE,col.names=FALSE)


run_t<-round(proc.time()-ptm,4) #Calculationg run time
print(paste("Runtime in seconds:",run_t[3]))
print(paste("Speed index (calc time for 500 permutations):",run_t[3]*500/arg$permutations/length(columns)))

print("Releasing cores")
stopCluster(cl) # Releasing acquired CPU cores
