setwd("c:/Users/Udi/SkyDrive/TCGA_CURATED/Rips")
arg<-list(10,NULL,NULL,NULL,NULL,"LUAD.h5","all",20,detectCores(),FALSE,TRUE,NULL,0.06,100,"syn","Annotations.csv",FALSE,FALSE,0,"PROCESSED_MAF_STADTRIM_2015-12-09.maf")
names(arg)<-c("epsilon","cut","topgenes","network","scan","matrix","columns","permutations","cores","log2","fdr","chunk","samples_threshold","g_score_threshold","score_type","anno","mutload","syn_control","rescale","maf")


#########g_scores#############
g_score_calc<-function(score_type,samples,genes,bin_matrix) { 
  mat_non_syn<-mat_non_syn[samples,genes,drop=FALSE]
  mat_syn<-mat_syn[samples,genes,drop=FALSE]
  mat_non_syn_bin<-bin_matrix[samples,genes,drop=FALSE]
  
    #G_scores type 1 - Based on  non-syn/(syn+non-syn) ratio
    syn_ratio<-colSums(mat_non_syn)/(colSums(mat_syn)+colSums(mat_non_syn)) #G_Score
    syn_ratio[syn_ratio=="NaN"]<-0 #Fixing 0 mutations columns 3
    g_score<-syn_ratio   

  
  return(g_score)
}


perm_values<-function(dict_matrix,column,matrix) {
  #Takes dictionary matrix and all samples- returns translated_matrix with corersponding values
  
  y<-apply(dict_matrix,2,function(perm_column) perm_column<-matrix[perm_column,column])
  
}



c_calc_fast2<-function(e_matrix)
  #Calculate connectivity value of a prticular column, 
  #c=sigma(i<>j) [ei*Aij*ej/ei*ej]
{
  c_base_outer_product<-e_matrix[,1] %o% e_matrix[,1] 
  c_base<-sum(c_base_outer_product[lower.tri(c_base_outer_product)]) # Calculating sum[e(i)*e(j)] for i!=j
    
  c_vector<-apply(e_matrix,2,function(e_column) {
    c<-sum(e_column[edges1]*e_column[edges2])
  })
  c_vector<-c_vector/c_base 
} 



connectivity_analysis<-function(columns_of_interest,bin_matrix,perm_dict) {
  
  print (paste("Preparing parallel environment. Acquiring",arg$cores,"Cores"))
  cl <- makeCluster(as.numeric(arg$cores))
  varlist=c("file_prefix","c_calc_fast2","edges1","edges2","samples","permutations","num_nodes","perm_values","arg","nodes","bin_matrix","columns","perm_dict")
  clusterExport(cl=cl, varlist=varlist,envir=environment())
  
  columns<-seq_along(columns_of_interest) #Subsetting columns
  split.column<-split(columns,ceiling(seq_along(columns)/arg$chunk))
  
  #count<-0
  #ans<-parLapply(cl,split.column,function (columns_range)  {
  ans<-lapply(split.column,function (columns_range)  {
    #calculating c-scores and p-values for each chunk of columns
    matrix2<-bin_matrix[,columns_range,drop=FALSE]
    
    perm_values_list<-as.list(rep(NA,length(columns_range))) #Creating list of "column" elements
    
    for (column in seq_along(columns_range)) {
      perm_values_list[[column]]<-perm_values(perm_dict,column,matrix2)
    }
      
    
    e_list<-perm_values_list
    #pi_list<-lapply(e_list,function(x) pii_matrix(as.matrix(x))) #columns_range elements in the list. Each element is a matrix representing pi_values of a gene. rows are nodes, columns are permutations. first column is non permuted.
    #c_vec_list<-lapply(pi_list,function (pi_matrix) c_calc_fast(as.matrix(pi_matrix)))
    c_vec_list<-lapply(e_list,function (e_matrix) c_calc_fast2(as.matrix(e_matrix)))
    
    #e_mean<-sapply(e_list,function (x) mean(x[,1])) #Taking mean of the first column (not permutations)
    #e_sd<-sapply(e_list,function (x) sd(x[,1]))
    #pi_values<-sapply(pi_list,function (x) x[,1]) #Is a matrix,each row is a node, each column in the matrix is pi values of a gene across nodes.
    #pi_frac<-apply(pi_values,2,function (x) sum(x!=0)/length(x))
    n_samples<-apply(matrix2,2,function (x) sum(x!=0))
    c_value<-sapply(c_vec_list,function (c_vec) c_vec[1]) #the first position is the connectivity value, the later are c_value for each permutation
   # p_value<-sapply(c_vec_list,function(c_vec) {
    #  p_value<-sum(c_vec>c_vec[1])/permutations})
    
    Genes<-colnames(matrix2)  
    results<-cbind(Genes,c_value,n_samples) #The variable names should match info_cols
    c_vec_df<-t(as.data.frame(c_vec_list))
    rownames(c_vec_df)<-Genes
    #write.csv(c_vec_df,paste0("output_",columns_range[1],".csv"))
    
    #Ans is columns_range length list. Each element contain to variables.
    #First variabl is "output" which is results matrix. the second element is pii_values matrix, its columns are genes and rows are nodes 
    
    ans<-list(results,c_vec_df)
    
    return(ans) 
    
  })
  print("Releasing cores")
  stopCluster(cl) # Releasing acquired CPU cores
  return(ans)
}


results_file<-function(ans) {
  if (length(columns_of_interest)==0) {
    final_results<-data.frame(Gene_Symbol="None_passed_threshold",p_value=NA,q_value=NA,q_value_integrated=NA)
    return(final_results)
  }
  final_results<-NULL
  for (i in 1:length(ans))
    final_results<-rbind(final_results,ans[[i]][[1]]) #Extracting output from ans
  
  final_results<-as.matrix(final_results,rownames.force = T,drop=FALSE)
  #pi_zero_genes<-final_results[,"pi_frac"]==0 #Probably not needed since all genes like that are out with g_score filtering
  #final_results<-final_results[!pi_zero_genes,,drop=FALSE]
  #final_results[pi_zero_genes,"p_value"]<-NA
  #q_value<-p.adjust(final_results[,"p_value"],"fdr")
  g_value<-columns_of_interest[rownames(final_results)]
  final_results<-cbind(final_results,g_value)
  #colnames(final_results)[grep("g_value",colnames(final_results))]<-paste0("g_score_",arg$score_type)
  #final_results<-final_results[order(final_results[,"q_value"]),,drop=FALSE]
  
  #if (arg$rescale==0) {final_results[,"Genes"]<-substring(final_results[,"Genes"],5)}
  if (arg$mutload==FALSE) {
    Gene_Symbol<-sapply(strsplit(final_results[,1],"|",fixed = TRUE),"[[",1)
    EntrezID<-sapply(strsplit(final_results[,1],"|",fixed = TRUE),"[[",2)
    final_results<-final_results[,-1,drop=FALSE] #Removing old genes column
    final_results<-cbind(Gene_Symbol,EntrezID,final_results)
  }
  return(final_results)
}


rescale<-function (cut,maf,mat_syn,mat_non_syn) {
  #This function input is cut value (number of samples above which to subsumple in log10 scale) , maf file , and the current mat_syn and mat-non_syn intersected matrices
  #The output is a subsumpled mat_non_syn_bin matrix
  cut<-10^cut #Reverting log10
  print ("Initiating rescaling process")
  mat_total<-mat_non_syn+mat_syn #Total number of point mutations before rescaling
  mutLoad<-rowSums(mat_total)
  above_cut<-mutLoad[mutLoad>cut]
  below_cut<-mutLoad[mutLoad<=cut]
  median(above_cut)
  median(below_cut)
  scale<-floor(median(above_cut)/median(below_cut)) # Rescaling parameter
  
  #Detecting mutloadmutated samples and removing from maf file based on scale
  maf<-read.delim(arg$maf,header = TRUE,as.is=T,comment.char = "#",sep="\t")
  maf$Tumor_Sample_Barcode<-substring(maf$Tumor_Sample_Barcode,1,15)
  
  mutloadmutated<-names(above_cut) #samples that appear to be hyper mutated
  x<-maf[maf$Tumor_Sample_Barcode %in% mutloadmutated,]
  shuffled_rows<-sample(1:nrow(x)) # Shuffling before subsampling
  x<-x[shuffled_rows,]
  rows_index_to_keep<-seq.int(1,nrow(x),by=(scale)) # Removing every row such that the ratio is 
  rows_index_to_remove<-setdiff(1:nrow(x),rows_index_to_keep)
  rownames_to_remove<-rownames(x[rows_index_to_remove,])
  maf<-maf[-match(rownames_to_remove,rownames(maf)),]
  maf$Column_name<-paste0("mut_",maf$Column_name) #Necessary for downstream process
  maf<-maf[maf$Tumor_Sample_Barcode %in% all_samples,] #Intersecting MAF file with all_samples
  
  #Creating rescaled matrices 
  all_genes<-sort(unique(maf$Column_name))
  all_samples<-sort(unique(maf$Tumor_Sample_Barcode))
  mat_syn<-matrix(0,length(all_samples),length(all_genes))  
  dimnames(mat_syn)<-list(all_samples,all_genes)
  mat_non_syn<-mat_syn #Replicating mat_syn
  
  #Creating table of Synonymous mutations for Sample vs Entrez_Gene_Id
  t_syn<-with(maf[maf$Synonymous,],table(Tumor_Sample_Barcode,Column_name))
  t_non_syn<-with(maf[!maf$Synonymous,],table(Tumor_Sample_Barcode,Column_name))
  
  #Plugging tables into matrices
  mat_syn[rownames(t_syn),colnames(t_syn)]<-t_syn
  mat_non_syn[rownames(t_non_syn),colnames(t_non_syn)]<-t_non_syn
  mat_syn<-mat_syn[sort(rownames(mat_syn)),sort(colnames(mat_syn))]
  mat_non_syn<-mat_non_syn[sort(rownames(mat_non_syn)),sort(colnames(mat_non_syn))]
  #Binary matrix for connectivity score
  mat_non_syn_bin_new<-ifelse(mat_non_syn>0,1,0) #Non synonymous binary matrix - will be used as input for c_score
  rescaled_matrices<-list(mat_non_syn_bin=mat_non_syn_bin_new,mat_syn=mat_syn,mat_non_syn=mat_non_syn)
  return(rescaled_matrices)
}




#rm(list=ls())

#del<-list.files(pattern = "*.csv")
#unlink(del)
#source('connectivity_functions2.R')
############################LOADING LIBRARIES############################

#setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/DATA/Agilent")
suppressWarnings({
  suppressMessages ({
    require(igraph,quietly = T,warn.conflicts = FALSE)
    library(rgexf,quietly = T,warn.conflicts = FALSE)
    library(jsonlite,quietly = T,warn.conflicts = FALSE)
    library(parallel,quietly = T,warn.conflicts = FALSE)
    library(getopt,quietly = T,warn.conflicts = FALSE)
    library(data.table,quietly = T,warn.conflicts = FALSE)
    library(rhdf5,quietly = T,warn.conflicts = FALSE)
    library(ggplot2,quietly = T,warn.conflicts = FALSE)
    library(stringr,quietly = T,warn.conflicts = FALSE)
    
    require("reshape2")
    require("bioDist")
  })
  
})


############################COMMAND LINE PARSING############################


#Setting defaults for debug mode
#arg<-list(20,NULL,NULL,NULL,NULL,"LUAD.h5","all",100,detectCores(),FALSE,TRUE,NULL,0.06,100,"syn","Annotations.csv",FALSE,FALSE,0,"PROCESSED_MAF_STADTRIM_2015-12-09.maf")
#names(arg)<-c("epsilon","cut","topgenes","network","scan","matrix","columns","permutations","cores","log2","fdr","chunk","samples_threshold","g_score_threshold","score_type","anno","mutload","syn_control","rescale","maf")

#Argument section handling
spec = matrix(c(
  "network", "n", 1, "character",
  "matrix", "m",1,"character",
  "columns", "i",1,"character",
  "g_score_threshold", "g",1,"integer",
  "permutations","p",2,"integer",
  "cores","q",1,"integer",
  "samples_threshold","t",1,"numeric", # Percentage of samples_of_interest
  "log2","l",2,"logical",
  "fdr","f",2,"logical",
  "chunk","k",2,"integer",  
  "score_type","s",1,"character",
  "anno","a",2,"character",
  "mutload","h",2,"logical",
  "syn_control","z",2,"logical",
  "rescale","r",2,"numeric",
  "maf","x",2,"character",
  "scan","y",2,"integer",
  "epsilon","e",2,"integer",
  "distance","d",2,"character",
  "topgenes","T",2,"integer",
  "cut","c",2,"numeric"
), byrow=TRUE, ncol=4)

#arg<-getopt(spec) #Conmment this line for debug mode


if ( is.null(arg$permutations ) ) {arg$permutations= 20}
if ( is.null(arg$log2 ) ) {arg$log2= FALSE}
if ( is.null(arg$fdr ) ) {arg$fdr= TRUE}
if ( is.null(arg$cores ) ) {arg$cores= 4}
if ( is.null(arg$chunk ) ) {arg$chunk= 25}
if ( is.null(arg$columns ) ) {arg$columns= "all"}
if ( is.null(arg$samples_threshold ) ) {arg$samples_threshold= 0.05}
if ( is.null(arg$anno ) ) {arg$anno= "Annotations.csv"}
if ( is.null(arg$mutload ) ) {arg$mutload= FALSE}
if ( is.null(arg$g_score_threshold ) ) {arg$g_score_threshold= 100}
if (arg$mutload==TRUE) {
  arg$samples_threshold<-0
  arg$g_score_threshold<-1
  arg$syn_control<-FALSE
}
if ( is.null(arg$syn_control ) ) {arg$syn_control= FALSE}
if ( is.null(arg$score_type ) ) {arg$score_type= "syn"}
if ( is.null(arg$rescale ) ) {arg$rescale= 0} else {
  if (is.null(arg$maf)) stop("PROCESSED MAF file must be provided for resscaling")
}

if ( is.null(arg$scan ) ) {arg$scan= 0}
if ( !is.null(arg$network ) ) {arg$scan= 0}

if ( is.null(arg$epsilon ) ) {arg$epsilon= 20}
if ( is.null(arg$topgenes ) ) {arg$topgenes= 2000}
if ( is.null(arg$cut) ) {arg$cut= 0.5}
if ( is.null(arg$dist) ) {arg$dist= 1}

PROJECT_NAME<-as.character(str_match(string = arg$matrix,pattern="\\w+.h5"))
PROJECT_NAME<-substring(text = PROJECT_NAME,first = 1,nchar(PROJECT_NAME)-3)
#Printing run parameters
print (paste("Number of permutations:",arg$permutations))
print (paste("Number of cores:",arg$cores))
print (paste("Chunk size:",arg$chunk))



####################################### MUTATIONS MATRIX HANDLING #############################


#Loading matrix file to memory and log transforming if log2=TRUE

#Loading matrix file to memory and log transforming if log2=TRUE
h5file<-arg$matrix
print (paste0("Loading ",h5file, " file to memory"))
mat_non_syn_bin<-h5read(h5file,"Mutations_Binary")
mat_non_syn<-h5read(h5file,"Mutations_NS")
mat_syn<-h5read(h5file,"Mutations_S")
mat_tpm<-h5read(h5file,"TPM")
all_samples<-h5read(h5file,"Mutations_Samples")
all_genes<-h5read(h5file,"Mutations_Genes")
tpm_genes<-h5read(h5file,"TPM_Genes")


rownames(mat_non_syn_bin)<-all_samples
colnames(mat_non_syn_bin)<-all_genes
rownames(mat_non_syn)<-all_samples
colnames(mat_non_syn)<-all_genes
rownames(mat_syn)<-all_samples
colnames(mat_syn)<-all_genes
rownames(mat_tpm)<-all_samples
colnames(mat_tpm)<-tpm_genes



guid<-round(runif(1, min = 300000, max = 399999),0)

##################RESCALING####################





if (arg$rescale!=0) { 
  subsampled_matrices<-rescale(arg$rescale,arg$maf,mat_syn,mat_non_syn)
  mat_non_syn_bin<-subsampled_matrices$mat_non_syn_bin
  mat_non_syn<-subsampled_matrices$mat_non_syn
  mat_syn<-subsampled_matrices$mat_syn
}

if (arg$log2==TRUE) {mat_tpm<-(2^mat_tpm)-1} #Preparing for calculation if matrix is log scale

#Info_cols is used to set information columns in results output file as well as names for the variables that constitutes those columns
#info_cols<-t(c("Genes","c_value","p_value","pi_frac","n_samples","e_mean","e_sd")) 






#################################################################
###############Generating mutational load histogra##############
#################################################################



mutload_matrix<-mat_non_syn+mat_syn #Total number of point mutations
mutload_dist<-rowSums(mutload_matrix)
if (arg$rescale!=0) {
  #png("hist_mutLoad_Rescaled.png")
  #hist(log10(mutload_dist),breaks = 100,main="After rescaling")
  #invisible(dev.off())  
  qplot(log10(mutload_dist),main = paste ("After rescale",PROJECT_NAME),) + ggsave(paste0("hist_mutLoad_Rescaled_",as.character(arg$rescale),"_",guid,".png"))
} else {
  #png("hist_mutLoad_NoRescaling.png")
  #hist(log10(mutload_dist),breaks = 100,main="Before rescaling")
  #invisible(dev.off())
  qplot(log10(mutload_dist),main = paste ("Before rescaling",PROJECT_NAME)) + ggsave(paste0("hist_mutLoad_NORescaling_",guid,".png"))
}




#################################################################################
####################### Functions Section #######################################
#################################################################################



#################################################################################################
################################ END OF FUNCTIONS SECTION########################################
##################################################################################################



som<-function(x) {sd(x)/mean(x)}
exp_mean<-colMeans(mat_tpm)
qplot(exp_mean,binwidth=0.1) + ggtitle("bin=0.1") + ggsave(paste0("exp_hist_",guid,".png"))

exp_var<-apply(mat_tpm[,exp_mean>arg$cut],2,som) #To avoid distortion, Calculate Coefficient of variationonly for genes with mean expression larger then arg$cut
topgenes<-names(head(sort(round(exp_var,3),decreasing = T),arg$topgenes))

exp_topgenes<-mat_tpm[,topgenes]
exp_topgenes<-t(exp_topgenes)
#cor_exp_topgenes<-1-cor(exp_topgenes)
if (arg$dist==1) {
	cor_exp_topgenes<-1-cor(exp_topgenes,method="pearson")
} 

if (arg$dist==2) {
	cor_exp_topgenes<-1-cor(exp_topgenes,method="spearman")
} 

if (arg$dist==3) {
	cor_exp_topgenes<-as.matrix(dist(t(exp_topgenes),upper = T))
}

#require("hopach")
#cor_exp_topgenes<-1-as.matrix(mutualInfo(t(exp_topgenes)))
  #as.matrix(distancematrix(t(exp_topgenes),"eucli"))
#g<-1-cor(exp_topgenes)
#sum(is.na(cor_exp_topgenes))



distance_set<-cor_exp_topgenes[lower.tri(cor_exp_topgenes,diag=FALSE)] #Recording lowe triangle values only.
epsilon_set<-sort(unique(distance_set)) #Extracting all possible distances
#EQUALIZE EPSILON - 
#epsilon_set<-seq(min(epsilon_set),max(epsilon_set),length.out = arg$epsilon)
epsilon_sub_index<-floor(seq(1,length(epsilon_set),length.out = arg$epsilon))
epsilon_set<-epsilon_set[epsilon_sub_index]
#epsilon_set<-round(epsilon_set,3)


#epsilon_min<-min(cor_exp_topgenes)
#epsilon_max<-max(cor_exp_topgenes)
#epsilon_set<-round(seq(epsilon_min,epsilon_max,length.out=arg$epsilon),2)


scan<-data.frame(networks=epsilon_set,resolution=NA,gain=NA,original_samples=NA,first_connected_samples=NA,edges_num=NA,samples_threshold=NA,above_samples_threshold=NA,above_gscore_threshold=NA,p_0.05=NA,q_0.1=NA,q_0.15=NA,q_0.2=NA,mutload=NA,parsing_time=NA,connectivity_time=NA,uid=NA,stringsAsFactors = F)
  
if (arg$scan!=0) { #Removing network files file for test mode
  scan<-scan[1:arg$scan,]
}





############################################################################################################
################### Running connectivity analysis for networks in scan table ##############################
############################################################################################################


dict<-matrix(NA,length(all_samples),arg$permutations+1) #rows= unique_sample_id cols= permutation ID, flash=permuted sample ID
#perm_values_list<-as.list(rep(NA,length(columns_of_interest))) #Creating list of "column" elements
perm_dict<-apply(dict,2,function(x) x<-sample(1:length(all_samples)))
perm_dict[,1]<-1:length(all_samples)


count<-0
c_matrix_list<-list()
#for (file in scan$networks[scan$mutload_connectivity]) {
for (epsilon in scan$networks) {
  count<-count+1  
  print("*********************************************")
  print (Sys.time())
  print (paste("Analyzing network:",epsilon,"-",count,"out of",nrow(scan)))
  print (paste("Number of permutations:",arg$permutations))
  print (paste("Number of CPU cores:",arg$cores))
  
  adj_mat<-ifelse(cor_exp_topgenes<=epsilon,1,0) 
  adj_mat[lower.tri(adj_mat,diag = TRUE)]<-0
  rownames(adj_mat)<-1:nrow(adj_mat)
  colnames(adj_mat)<-1:ncol(adj_mat)
  graph_igraph<-graph.adjacency(adj_mat,diag = F) #Converting the adjacency matrix into igraph format
  edges<-get.edgelist(graph_igraph,names=FALSE) #Extracting edges from the graph
  nodes<-lapply(1:nrow(adj_mat),function (x) x) #Generating node list - every node is a sample
  names(nodes)<-1:nrow(adj_mat)
  
  edges_num<-nrow(edges) #Recording number of edges
  print (paste("Number of edges in network:",edges_num))
  scan[scan$networks==epsilon,]$edges_num<-edges_num
  
  samples<-unique(unlist(nodes)) #Recording samples in all nodes - here every node==sample
  matrix1<-mat_non_syn_bin[samples,,drop=FALSE] #Subseting matrix to contain only samples in first connected graph 
  
  #Extracting columns from arguments
  columns<-seq_along(colnames(matrix1))
  #Removing columns below samples_threshold from the first connected graph
  #matrix1<-matrix1[,columns,drop=FALSE] #Subsetting for selected columns
  
  #samples_of_interest<-rownames(matrix1)
  samples_of_interest<-rownames(matrix1)
  
  # Taking record of sample sizes
  
  scan[scan$networks==epsilon,]$original_samples<-length(all_samples)
  scan[scan$networks==epsilon,]$first_connected_samples<-length(samples_of_interest)
  
  
  
  ##########################################
  ############THRESHOLDING SECTION###########
  ##########################################
  
  
  # selecting genes based on thresholds
  
  samples_threshold<-ceiling(arg$samples_threshold*length(samples_of_interest))
  scan[scan$networks==epsilon,]$samples_threshold<-samples_threshold
  print (paste("Samples in original dataset:",length(all_samples)))
  print (paste("Samples in first connected graph:",length(samples_of_interest)))
  print (paste("Samples threshold is set to:",samples_threshold))
  print (paste("Top genes score threshold is set to:",arg$g_score_threshold))
  
  genes_number_of_samples<-apply(matrix1,2,function (x) sum(x!=0)) #Counting non_zero samples for each column
  #genes_below_samples_threshold<-names(which(genes_number_of_samples<samples_threshold))
  genes_above_samples_threshold<-names(which(genes_number_of_samples>=samples_threshold)) #For filtering by number  of mutations exist in a sample
  g_score<-g_score_calc(arg$score_type,samples_of_interest,genes_above_samples_threshold,matrix1) #Syn/old only over sample thresholded genes
  
  columns_of_interest<-head(sort(g_score,decreasing = T),arg$g_score_threshold) #Filtering by g-score
  
  print(paste0("Columns above threshold: ",length(columns_of_interest)))
  
  scan[scan$networks==epsilon,]$above_samples_threshold<-length(genes_above_samples_threshold)
  scan[scan$networks==epsilon,]$above_gscore_threshold<-length(columns_of_interest)
  
  
  matrix1<-matrix1[,names(columns_of_interest),drop=FALSE] #Subsetting matrix to have above threshold columns
  
  if (arg$mutload==TRUE) {
    #Adding to matrix1 a column with mutation rate, this will be used to assess mutload mutated samples. 
    
    matrix1<-as.matrix(mutload_dist[samples_of_interest],drop=FALSE)
    colnames(matrix1)<-"mutLoad"
    columns_of_interest<-"mutLoad"
  }
  
  #Initializing results file name and unique id
  unique_id<-round(runif(1, min = 111111, max = 222222),0)
  file_prefix<-paste0(epsilon,"_",PROJECT_NAME,"-",unique_id,"_",guid,"-",Sys.Date())
  print(paste("File unique identifier:",unique_id))
  scan[scan$networks==epsilon,]$uid<-unique_id
  
  
  
    
  #Writing log file
  suppressWarnings(write.table(as.character(arg) ,paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("Number of permutations: ",arg$permutations),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("Samples threshold: ",arg$samples_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("g_score threshold: ",arg$g_score_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("Columns above threshold:",length(columns_of_interest)),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste0("Original sample size:",length(all_samples)),paste0(file_prefix,"_log.csv"),append=TRUE))
  #suppressWarnings(write.table(paste0("First connected sample size:",length(samples_of_interest)),paste0(file_prefix,"_log.csv"),append=TRUE))
  
  
  
  permutations<-arg$permutations
  edges1<-edges[,1] #Nodes i
  edges2<-edges[,2] #Nodes j
  num_nodes<-length(nodes)
  
  
  
  ptm<-proc.time()
  
  
  
  
  ######################### CONNECTIVITY ANALYSIS ##############################################3
  
  #if (edges_num!=0) 
    
    
    print ("Starting connectivity analysis:")
    if (arg$mutload==TRUE) {
      print (paste("Mutload Connectivity for Graph"))
      ans<-connectivity_analysis(columns_of_interest,matrix1,perm_dict)
     
      
      final_results<-results_file(ans)
      scan[scan$networks==epsilon,]$mutload<-final_results[,"p_value"]
      
      
    } else {   # Genes analysys
      print (paste("Genes Connectivity for Graph"))
      ans<-connectivity_analysis(columns_of_interest,matrix1,perm_dict)
      final_results<-results_file(ans)
      #scan[scan$networks==epsilon,]$p_0.05<-sum(final_results[,"p_value"]<=0.05)
      #scan[scan$networks==epsilon,]$q_0.1<-sum(final_results[,"q_value"]<=0.1)
      #scan[scan$networks==epsilon,]$q_0.15<-sum(final_results[,"q_value"]<=0.15)
      #scan[scan$networks==epsilon,]$q_0.2<-sum(final_results[,"q_value"]<=0.2)
    }
    
    
    c_matrix<-NULL #rows:genes, cols: permutations, value: connectivity value.
    for (i in 1:length(ans)) {
      c_matrix<-rbind(c_matrix,ans[[i]][[2]])
    }
    c_matrix_file<-paste0(file_prefix,"_c_matrix.csv")
    write.csv(c_matrix,c_matrix_file)
    
    ################ pii_value generations goes here - pushed to end##################
    c_matrix_list[[as.character(epsilon)]]<-c_matrix     #Recordind c_matrix
    
    
    
    ##################Synonymous control###############
    if (arg$syn_control==TRUE & length(columns_of_interest)!=0) {
      print ("Starting control connectivity analysis:")
      matrix1<-ifelse(mat_syn>0,1,0) #Using synonymous matrix as reference
      matrix1<-matrix1[samples_of_interest,names(columns_of_interest),drop=FALSE] # Subsetting for samples of interest
      ans<-connectivity_analysis(columns_of_interest,matrix1,perm_dict) #Running connectivity analysis
      final_results_control<-results_file(ans)
      
      #Coercing non_syn and control results
      final_results_control<-final_results_control[,c("n_samples","c_value","p_value","q_value"),drop=FALSE]
      colnames(final_results_control)<-c("n_samples_con","c_value_con","p_value_con","q_value_con")
      missing_genes<-setdiff(rownames(final_results),rownames(final_results_control)) #Genes that do not exist in final_Results needed to be completed with NA and zeros
      
      missing_n_samples_con<-colSums(matrix1[,missing_genes])
      missing_p_value_con<-rep(NA,length(missing_genes))
      missing_q_value_con<-rep(NA,length(missing_genes))
      missing_c_value_con<-rep(NA,length(missing_genes))
      missing_x<-data.frame(missing_n_samples_con,missing_c_value_con,missing_p_value_con,missing_q_value_con)
      final_results_control<-rbind(final_results_control,as.matrix(missing_x))
      final_results_control<-final_results_control[rownames(final_results),,drop=FALSE]
      final_results<-cbind(final_results,final_results_control)
      final_results<-final_results[,c("Gene_Symbol","EntrezID","c_value","p_value","n_samples","q_value","c_value_con","p_value_con","n_samples_con","g_score_syn"),drop=FALSE]
      
      #Calculating integrated p_value
      n<-as.numeric(final_results[,"n_samples"])
      n_con<-as.numeric(final_results[,"n_samples_con"])
      p<-as.numeric(final_results[,"p_value"])
      p_con<-as.numeric(final_results[,"p_value_con"])
      
      p_integrated<-p_integrate(p,p_con,n,n_con)
      q_integrated<-p.adjust(p_integrated,"fdr")
      
      final_results<-cbind(final_results,p_integrated,q_integrated)
      
    }
    
    ##################Fixing final results file and some printings#####################
    
    
    
    if (arg$mutload==TRUE) {
      file_sufix<-"_mutload_results.csv"
    } else {file_sufix<-"_genes_results.csv"}
    
    if (arg$mutload==FALSE) {
      final_results[,"Gene_Symbol"]<-sapply(final_results[,"Gene_Symbol"],function (x) strsplit(x,"mut_")[[1]][2])  
    }
    
    write.table(final_results,paste0(file_prefix,file_sufix),row.names=FALSE,sep=",")
    
    run_t<-round(proc.time()-ptm,4) #Calculating run time
    scan[scan$networks==epsilon,]$connectivity_time<-run_t[3]
    speed_index<-run_t[3]*500/arg$permutations/length(columns)
    print(paste("Runtime in seconds:",run_t[3]))
    print(paste("Speed index (calc time for 500 permutations):",speed_index)) 
    
  #}
   # else {
    #print ("NO EDGES IN THE GRAPH ONLY SCAN FILE IS PRODUCED")
    #scan[scan$networks==epsilon,][,7:ncol(scan)]<-0  # NO EDGES IN THE GRAPH
    
  #}
  
  
  
}



#Writing scanner summary file:
write.csv(scan,paste0("scan_summary_",guid,".csv"))



###############MOVING FILES TO Results DIR####################



if (arg$mutload==TRUE) {
  results_dir<-paste0("results_",PROJECT_NAME,"_",guid,"_mutload")
} else { results_dir<-paste0("results_",PROJECT_NAME,"_",guid,"_genes")}

if (arg$rescale!=0) {results_dir<-paste0(results_dir,"_rescaled_",arg$rescale)}

results_tar<-paste0(results_dir,".tar.gz")  


print (paste("Moving files to Results Directory:",results_dir))

dir.create(results_dir)
#if (file.exists("Rplots.pdf")) {file.remove("Rplots.pdf")}

files_to_move_to_results<-list.files(pattern = as.character(guid))
x<-file.rename(files_to_move_to_results,paste0(results_dir,"/",files_to_move_to_results))
tar(results_tar,results_dir,compression="gzip")



x<-list()
for (i in seq_along(columns_of_interest)) {
  gene<-names(columns_of_interest[i])
  x[[gene]]<-sapply(c_matrix_list,function(epsilon) epsilon[gene,])
}




for (i in 2:nrow(x[[36]])) {
  plot(epsilon_set,x[[36]][i,],type="l",xlab="Epsilon",ylab=("Connectivity"))
  par(new=TRUE)
}
plot(epsilon_set,x[[36]][1,],col="red",type="l",lwd=10,xlab="Epsilon",ylab=("Connectivity"))


##qplot(epsilon_set,x[[2]][1,],geom="line",color="white",size=7)
#qplot(epsilon_set,x[[1]][1,]) + geom_line(color="blue",size=5)

dist_mean<-function(c_dist) {
  average_slope<-c_dist[1]*epsilon_set[1]
  for (i in (2:length(c_dist))) {
    dslope<-(c_dist[i]-c_dist[i-1])*epsilon_set[i]
    average_slope<-average_slope+dslope
  }
  return(average_slope)
}

dist_mean2<-function(c_dist) {
  average_slope<-c_dist[1]*1
  for (i in (2:length(c_dist))) {
    dslope<-(c_dist[i]-c_dist[i-1])*i
    average_slope<-average_slope+dslope
  }
  return(average_slope)
}

#dist_area<-function(c_dist) {
#  a<-0
#  for (i in (2:length(c_dist))) {
#    rect_area<-(epsilon_set[i]-epsilon_set[i-1])*(c_dist[i-1]+c_dist[i])/2
#    a<-a+rect_area
#  }
#  return(a)
#}


dist_mean_list<-lapply(x,function (z) apply(z,1,dist_mean))
#dist_mean_list<-lapply(x,function (z) apply(z,1,dist_area))


p_value<-sapply(dist_mean_list,function(q) {
  p_value<-sum(q<q[1])/permutations})

q_value<-p.adjust(p_value,"fdr")
genes<-sort(names(columns_of_interest))


results<-data.frame(p_value[genes],q_value[genes])
results<-results[order(results$q_value.genes.),]
head(results,20)

write.csv(results,paste0(file_prefix,"_new_results.csv"))

