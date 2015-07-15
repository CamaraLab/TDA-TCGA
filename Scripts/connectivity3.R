#setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/DATA")
#Preping environment, loading necessary libraries

library(igraph)
library(rgexf)
library(jsonlite)
library(parallel)
library(getopt)
library(data.table)
library(rhdf5)

#Setting defaults for debug mode
arg<-list("LUAD_Cor_MDS_22_3_intersect","LUAD.h5","all",500,detectCores(),FALSE,TRUE,50,50,50,"lam")
names(arg)<-c("name","matrix","columns","permutations","cores","log2","fdr","chunk","samples_threshold","g_score_threshold","score_type")

#Argument section handling
spec = matrix(c(
  "name", "n", 1, "character",
  "matrix", "m",1,"character",
  "columns", "c",1,"character",
  "g_score_threshold", "g",1,"integer",
  "permutations","p",2,"integer",
  "cores","q",1,"integer",
  "samples_threshold","t",1,"integer", # Minimum number of samples to express column - Removing columns below that
  "log2","l",2,"logical",
  "fdr","f",2,"logical",
  "chunk","k",2,"integer",  
  "score_type","s",1,"character"
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode

if ( is.null(arg$permutations ) ) {arg$permutations= 500}
if ( is.null(arg$log2 ) ) {arg$log2= TRUE}
if ( is.null(arg$fdr ) ) {arg$fdr= TRUE}
if ( is.null(arg$cores ) ) {arg$cores= detectCores()}
if ( is.null(arg$chunk ) ) {arg$chunk= 200}
if ( is.null(arg$columns ) ) {arg$columns= "all"}
if ( is.null(arg$samples_threshold ) ) {arg$samples_threshold= 0}



#Loading matrix file to memory and log transforming if log2=TRUE
matrix_full_name<-arg$matrix
print (paste0("Loading ",matrix_full_name, " file to memory"))
matrix1<-h5read("LUAD.h5","Mutations_Binary")
mat_non_syn<-h5read("LUAD.h5","Mutations_NS")
mat_syn<-h5read("LUAD.h5","Mutations_S")
all_samples<-h5read("LUAD.h5","Mutations_Samples")
all_genes<-h5read("LUAD.h5","Mutations_Genes")
rownames(matrix1)<-all_samples
colnames(matrix1)<-all_genes
rownames(mat_non_syn)<-all_samples
colnames(mat_non_syn)<-all_genes
rownames(mat_syn)<-all_samples
colnames(mat_syn)<-all_genes
if (arg$log2==TRUE) {matrix1<-(2^matrix1)-1} #Preparing for calculation if matrix is log scale



column_range<-function(col_range)
  #Gets column range from arg$column and parse it. Also removes columns below threshold
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


#########g_scores#############
g_score<-function(score,samples,genes) {
  mat_non_syn<-mat_non_syn[samples,genes]
  mat_syn<-mat_syn[samples,genes]
  mat_non_syn_bin<-matrix1[samples,genes]
  if (score=="syn") {
    #G_scores type 1 - Based on  non-syn/(syn+non-syn) ratio
    syn_ratio<-colSums(mat_non_syn)/(colSums(mat_syn)+colSums(mat_non_syn)) #G_Score
    syn_ratio[syn_ratio=="NaN"]<-0 #Fixing 0 mutations columns 3
    g_score<-syn_ratio   
  }
  
  if (score=="lam") {
    # G_Scores type 2 - Based on gene lengths
    anno<-read.csv("Annotations.csv")
    Lg<-as.numeric(anno$length[match(names(genes),anno$Symbol)])
    names(Lg)<-names(genes)
    genes_with_known_length<-names(Lg[!is.na(Lg)])
    Lg<-Lg[genes_with_known_length] #Removing unknown length EntrezId's
    L<-sum(Lg) #Total Coding region length
    ns<-rowSums(mat_non_syn[,genes_with_known_length]) #Sum of non syn mutations for each row
    S_lambda<-rowSums(sapply(samples,function (x) {
      lambda<-Lg*ns[x]/L
      if(sum(lambda==0)) {ans<-lambda} else {
        ans<-(-1)*mat_non_syn_bin[x,genes_with_known_length]*log(1-exp(-lambda))  
      }
      return(ans)
    }))
    #suppressWarnings(S_lambda<-as.numeric(c(S_lambda,setdiff(colnames(mat_non_syn),names(S_lambda)))))
    g_score<-S_lambda
  }
  
  if (score=="old") {
    #G_Score 3
    #Divides each element by the sum of the corresponding row sum.
    # Returns zero in case of division by zero
    NS<-rowSums(mat_non_syn_bin)
    mat_NS<-mat_non_syn_bin/NS
    mat_NS[mat_NS=="NaN"]<-0
    NS<-colSums(mat_NS)
    g_score<-NS
  }
  H5close()
  return(g_score)
}


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

#Extracting columns from arguments
columns<-column_range(arg$columns)
#Removing columns below samples_threshold from the first connected graph
matrix1<-matrix1[,columns] #Subsetting for selected columns
genes_number_of_samples<-apply(matrix1,2,function (x) sum(x!=0)) #Counting non_zero samples for each column

genes_below_samples_threshold<-which(genes_number_of_samples<arg$samples_threshold)

genes_above_samples_threshold<-which(genes_number_of_samples>=arg$samples_threshold) #For filtering by number  of mutations exist in a sample

#Choosing genes based on score
samples_of_interest<-rownames(matrix1)

#print(head(sort(g_score(arg$score_type,samples_of_interest,genes_above_samples_threshold)),arg$g_score_threshold))
columns_of_interest<-head(sort(g_score(arg$score_type,samples_of_interest,genes_above_samples_threshold),decreasing=TRUE),arg$g_score_threshold) #Filtering by g-score
#columns_of_interest<-head(sort(g_score(1),decreasing=TRUE),100) #Filtering by g-score
#columns_of_interest<-head(sort(g_score(matrix1),decreasing=TRUE),arg$samples_threshold) #Filtering by g-score
print(paste0("Columns above threshold: ",length(columns_of_interest)))


matrix1<-matrix1[,names(columns_of_interest)] #Subsetting matrix to have above threshold columns
columns<-seq_along(columns_of_interest) #Subsetting columns



#Initializing rolling results file
unique_id<-round(runif(1, min = 111111, max = 222222),0)
file_prefix<-paste0(arg$name,"_",arg$matrix,"-",unique_id,"-",Sys.Date())
print(paste("File unique identifier:",unique_id))

#Info_cols is used to set inforation columns in output file as well as names for the variables that constitutes those columns

info_cols<-t(c("Genes","c_scores","p_values","pi_frac","n_samples","e_mean","e_sd")) 
write.table(info_cols,paste0(file_prefix,"_results_rolling.csv"),sep=",",col.names=FALSE,row.names=FALSE)

write.csv(names(genes_below_samples_threshold),paste0(file_prefix,"_thresholded_genes1.csv"))
write.csv(setdiff(names(genes_above_samples_threshold),names(columns_of_interest)),paste0(file_prefix,"_thresholded_genes2.csv"))

#Writing log file
suppressWarnings(write.table(as.character(arg) ,paste0(file_prefix,"_log.csv"),append=TRUE))
suppressWarnings(write.table(paste("Number of permutations: ",arg$permutations),paste0(file_prefix,"_log.csv"),append=TRUE))
suppressWarnings(write.table(paste("Samples threshold: ",arg$samples_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
suppressWarnings(write.table(paste("g_score threshold: ",arg$g_score_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
suppressWarnings(write.table(paste("Columns above threshold:",length(columns_of_interest)),paste0(file_prefix,"_log.csv"),append=TRUE))
#perm_values<-function(dict_matrix,samples_relabling_table,column,matrix) {
perm_values<-function(dict_matrix,column,matrix) {
  #Takes dictionary matrix and all samples- returns translated_matrix with corersponding values
  
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

split.column<-split(columns,ceiling(seq_along(columns)/arg$chunk))


print (paste("Preparing parallel environment. Acquiring",arg$cores,"Cores"))
cl <- makeCluster(as.numeric(arg$cores))
varlist=c("file_prefix","c_calc_fast","c_calc_fast","pii_matrix","e_matrix","edges1","edges2","samples","permutations","num_nodes","perm_values","arg","nodes","matrix1","largest_cluster_nodes","info_cols","columns","samples_relabling_table")
clusterExport(cl=cl, varlist=varlist,envir=environment())

print (Sys.time())
ptm<-proc.time()

ans<-parLapply(cl,split.column,function (columns_range)  {
  #calculating c-scores and p-values for each chunk of columns
  
  matrix2<-matrix1[,columns_range]
  dict<-matrix(NA,length(samples),permutations+1) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
  perm_dict_list<-as.list(rep(NA,length(columns_range))) #Creating list of "column" elements
  perm_values_list<-as.list(rep(NA,length(columns_range))) #Creating list of "column" elements
  
  perm_dict_list<-lapply(perm_dict_list,function(x) { #Each element in perm_dict list will get permutation matrix 
    x<-apply(dict,2,function(x) x<-sample(samples))
    x[,1]<-1:length(samples)
    return(x)
  }) 
  
  for (column in seq_along(columns_range)) {
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
  ans<-list(output,rbind(Genes,pi_values))
  return(ans) 
  
})


#Generating results table. 
final_results<-NULL
for (i in 1:length(ans))
  final_results<-rbind(final_results,ans[[i]][[1]]) #Extracting output from ans


pi_zero_genes<-final_results[,"pi_frac"]==0 #Probably not needed since all genes like that are out with g_Score filtering
final_results<-final_results[!pi_zero_genes,]
q_value<-p.adjust(final_results[,"p_values"],"fdr")
g_value<-columns_of_interest[rownames(final_results)]
final_results<-cbind(final_results,q_value,g_value)
colnames(final_results)[9]<-paste0("g_score_",arg$score_type)
final_results<-final_results[order(final_results[,"q_value"]),]



#Generating pii_values table
pi_values_table<-NULL
for (i in 1:length(ans))
  pi_values_table<-cbind(pi_values_table,ans[[i]][[2]]) #Extracting pii_values from ans
colnames(pi_values_table)<-pi_values_table[1,]
pi_values_table<-pi_values_table[,!pi_zero_genes]


#Writing final results and pii_values files
write.table(final_results,paste0(file_prefix,"_results_final.csv"),row.names=FALSE,sep=",")
write.table(pi_values_table,paste0(file_prefix,"_pii_values.csv"),sep=",",row.names=FALSE,col.names=FALSE)


run_t<-round(proc.time()-ptm,4) #Calculationg run time
speed_index<-run_t[3]*500/arg$permutations/length(columns)
print(paste("Runtime in seconds:",run_t[3]))

print(paste("Speed index (calc time for 500 permutations):",speed_index))

suppressWarnings(write.table(paste("Runtime in seconds:",run_t[3]),paste0(file_prefix,"_log.csv"),append=TRUE))
suppressWarnings(write.table(paste("Speed index:",speed_index),paste0(file_prefix,"_log.csv"),append=TRUE))

print("Releasing cores")
stopCluster(cl) # Releasing acquired CPU cores
