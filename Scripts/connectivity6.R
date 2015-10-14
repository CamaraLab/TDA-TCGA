
############################LOADING LIBRARIES############################

#setwd("c:/users/udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/DATA/Agilent")
#Preping environment, loading necessary libraries
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
  })
  
})


############################COMMAND LINE PARSING############################


#Setting defaults for debug mode
arg<-list("LUAD_Cor_PCA_15_1.5","LUAD.h5","all",200,detectCores(),FALSE,TRUE,NULL,1,100,"syn","Annotations.csv",FALSE,TRUE,0,"PROCESSED_hgsc.bcm.edu_COAD.IlluminaGA_DNASeq.1.somatic.v.2.1.5.0.maf")
names(arg)<-c("network","matrix","columns","permutations","cores","log2","fdr","chunk","samples_threshold","g_score_threshold","score_type","anno","mutload","syn_control","rescale","maf")

#Argument section handling
spec = matrix(c(
  "network", "n", 1, "character",
  "matrix", "m",1,"character",
  "columns", "c",1,"character",
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
  "scan","y",2,"integer"
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode

if ( is.null(arg$permutations ) ) {arg$permutations= 500}
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

if ( is.null(arg$syn_control ) ) {arg$syn_control= TRUE}
if ( is.null(arg$score_type ) ) {arg$score_type= "syn"}
if ( is.null(arg$rescale ) ) {arg$rescale= 0} else {
  if (is.null(arg$maf)) stop("PROCESSED MAF file must be provided for resscaling")
}

if ( is.null(arg$scan ) ) {arg$scan= 0}
if ( !is.null(arg$network ) ) {arg$test_mode= 0}


#Printing run parameters
print (paste("Number of permutations:",arg$permutations))
print (paste("Number of cores:",arg$cores))
print (paste("Chunk size:",arg$chunk))
      


####################################### MUTATIONS MATRIX HANDLING #############################


#Loading matrix file to memory and log transforming if log2=TRUE
h5file<-arg$matrix
print (paste0("Loading ",h5file, " file to memory"))
mat_non_syn_bin<-h5read(h5file,"Mutations_Binary")
mat_non_syn<-h5read(h5file,"Mutations_NS")
mat_syn<-h5read(h5file,"Mutations_S")
all_samples<-h5read(h5file,"Mutations_Samples")
all_genes<-h5read(h5file,"Mutations_Genes")


rownames(mat_non_syn_bin)<-all_samples
colnames(mat_non_syn_bin)<-all_genes
rownames(mat_non_syn)<-all_samples
colnames(mat_non_syn)<-all_genes
rownames(mat_syn)<-all_samples
colnames(mat_syn)<-all_genes




	##################RESCALING####################


	 if (arg$rescale!=0) {
	 cut<-10^arg$rescale
	 print ("Initiating rescaling process")
	 #maf<-read.delim("../../COAD_TEST/Mutations/PROCESSED_hgsc.bcm.edu_COAD.IlluminaGA_DNASeq.1.somatic.v.2.1.5.0.maf",header = TRUE,as.is=T,comment.char = "#",sep="\t")
	 maf<-read.delim(arg$maf,header = TRUE,as.is=T,comment.char = "#",sep="\t")
	 maf$Tumor_Sample_Barcode<-substring(maf$Tumor_Sample_Barcode,1,15)
	 mat_total<-mat_non_syn+mat_syn #Total number of point mutations
	 mutLoad<-rowSums(mat_total)
	 above_cut<-mutLoad[mutLoad>cut]
	 below_cut<-mutLoad[mutLoad<=cut]
	 median(above_cut)
	 median(below_cut)
	 scale<-floor(median(above_cut)/median(below_cut)) # Rescaling parameter

	 #Detecting mutloadmutated samples and removing from maf file based on scale
	 mutloadmutated<-names(above_cut)
	 x<-maf[maf$Tumor_Sample_Barcode %in% mutloadmutated,]
	 rows_to_keep<-rownames(x[seq.int(1,nrow(x),by =round(scale)),]) # Removing every SCALEth row
	 rows_to_remove<-setdiff(rownames(x),rows_to_keep)
	 maf<-maf[-match(rows_to_remove,rownames(maf)),]

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
	 mat_non_syn_bin<-ifelse(mat_non_syn>0,1,0) #Non synonymous binary matrix - will be used as input for c_score

	 }


if (arg$log2==TRUE) {mat_non_syn_bin<-(2^mat_non_syn_bin)-1} #Preparing for calculation if matrix is log scale

#Info_cols is used to set information columns in results output file as well as names for the variables that constitutes those columns
info_cols<-t(c("Genes","c_value","p_value","pi_frac","n_samples","e_mean","e_sd")) 






	#################################################################
	###############Generating mutational load histogra##############
	#################################################################


	 if (arg$mutload==TRUE) {
	 
     
		 
		  mutload_matrix<-mat_non_syn+mat_syn #Total number of point mutations
		  mutload_dist<-rowSums(mutload_matrix)
		  if (arg$rescale!=0) {
		   png("hist_mutLoad_Rescaled.png")
		   hist(log10(mutload_dist),breaks = 100,main="After rescaling")
		   invisible(dev.off())  
		  } else {
		   png("hist_mutLoad_NoRescaling.png")
		   hist(log10(mutload_dist),breaks = 100,main="Before rescaling")
		   invisible(dev.off())
		  }
		  
	 }




#################################################################################
####################### Functions Section #######################################
#################################################################################


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
g_score_calc<-function(score_type,samples,genes) { 
  mat_non_syn<-mat_non_syn[samples,genes,drop=FALSE]
  mat_syn<-mat_syn[samples,genes,drop=FALSE]
  mat_non_syn_bin<-matrix1[samples,genes,drop=FALSE]
  if (score_type=="syn") {
    #G_scores type 1 - Based on  non-syn/(syn+non-syn) ratio
    syn_ratio<-colSums(mat_non_syn)/(colSums(mat_syn)+colSums(mat_non_syn)) #G_Score
    syn_ratio[syn_ratio=="NaN"]<-0 #Fixing 0 mutations columns 3
    g_score<-syn_ratio   
  }
  
  if (score_type=="lam") {
    # G_Scores type 2 - Based on gene lengths
    anno<-read.csv(arg$anno)
    Lg<-as.numeric(anno$length[match(substring(genes,5),paste0(anno$Symbol,"|",anno$EntrezID))])
    names(Lg)<-genes
    genes_with_known_length<-names(Lg[!is.na(Lg)])
    Lg<-Lg[genes_with_known_length] #Removing unknown length EntrezId's
    L<-sum(Lg) #Total Coding region length
    ns<-rowSums(mat_non_syn[,genes_with_known_length]) #Sum of non syn mutations for each row
    S_lambda<-rowSums(sapply(samples,function (x) {
      lambda<-Lg*ns[x]/L
      ans<-(-1)*mat_non_syn_bin[x,genes_with_known_length]*log(1-exp(-lambda))  
      return(ans)
    }))
    g_score<-S_lambda
  }
  
  if (score_type=="old") {
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




connectivity_analysis<-function(columns_of_interest,matrix) {
  
  print (paste("Preparing parallel environment. Acquiring",arg$cores,"Cores"))
  cl <- makeCluster(as.numeric(arg$cores))
  varlist=c("file_prefix","c_calc_fast","c_calc_fast","pii_matrix","e_matrix","edges1","edges2","samples","permutations","num_nodes","perm_values","arg","nodes","matrix","largest_cluster_nodes","info_cols","columns","samples_relabling_table")
  clusterExport(cl=cl, varlist=varlist,envir=environment())
  
  columns<-seq_along(columns_of_interest) #Subsetting columns
  split.column<-split(columns,ceiling(seq_along(columns)/arg$chunk))
  
  
  ans<-parLapply(cl,split.column,function (columns_range)  {
    #ans<-lapply(split.column,function (columns_range)  {
    #calculating c-scores and p-values for each chunk of columns
    
    matrix2<-matrix[,columns_range,drop=FALSE]
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
    
    e_mean<-sapply(e_list,function (x) mean(x[,1])) #Taking mean of the first column (not permutations)
    e_sd<-sapply(e_list,function (x) sd(x[,1]))
    pi_values<-sapply(pi_list,function (x) x[,1]) #Is a matrix,each row is a node, each column in the matrix is pi values of a gene across nodes.
    pi_frac<-apply(pi_values,2,function (x) sum(x!=0)/length(x))
    n_samples<-apply(matrix2,2,function (x) sum(x!=0))
    c_value<-sapply(c_vec_list,function (c_vec) c_vec[1])
    p_value<-sapply(c_vec_list,function(c_vec) {
      p_value<-sum(c_vec>c_vec[1])/permutations})
    Genes<-colnames(matrix2)  
    
    output<-cbind(Genes,c_value,p_value,pi_frac,n_samples,e_mean,e_sd) #The variable names should match info_cols
    #write.table(output,paste0(file_prefix,"_results_rolling.csv"),append=TRUE,sep=",",col.names=FALSE,row.names=FALSE)
    
    #Ans is columns_range length list. Each element contain to variables.
    #First variabl is "output" which is results matrix. the second element is pii_values matrix, its columns are genes and rows are nodes 
    ans<-list(output,rbind(Genes,pi_values))
    
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
  pi_zero_genes<-final_results[,"pi_frac"]==0 #Probably not needed since all genes like that are out with g_score filtering
  #final_results<-final_results[!pi_zero_genes,,drop=FALSE]
  final_results[pi_zero_genes,"p_value"]<-NA
  q_value<-p.adjust(final_results[,"p_value"],"fdr")
  g_value<-columns_of_interest[rownames(final_results)]
  final_results<-cbind(final_results,q_value,g_value)
  colnames(final_results)[grep("g_value",colnames(final_results))]<-paste0("g_score_",arg$score_type)
  final_results<-final_results[order(final_results[,"q_value"]),,drop=FALSE]
  
  #if (arg$rescale==0) {final_results[,"Genes"]<-substring(final_results[,"Genes"],5)}
  if (arg$mutload==FALSE) {
    Gene_Symbol<-sapply(strsplit(final_results[,1],"|",fixed = TRUE),"[[",1)
    EntrezID<-sapply(strsplit(final_results[,1],"|",fixed = TRUE),"[[",2)
    final_results<-final_results[,-1,drop=FALSE] #Removing old genes column
    final_results<-cbind(Gene_Symbol,EntrezID,final_results)
  }
  return(final_results)
}


#P_values integration
p_integrate <- function (p,p_con,n,n_con)
  # Gets p_values together with popultaion size. Return integrated p_value
{
  z<-qnorm(p,lower.tail=FALSE)
  z_con<-qnorm(1-p_con,lower.tail=FALSE)
  z_weighted<-(n*z+pmin(n,n_con)*z_con)/sqrt(n^2+pmin(n,n_con)^2)
  p_weighted<-pnorm(z_weighted,lower.tail = FALSE)
  p_weighted[is.na(p_weighted) & !is.nan(p_weighted)]<-p[is.na(p_weighted) & !is.nan(p_weighted)]
  return (p_weighted)
}



#################################################################################################
################################ END OF FUNCTIONS SECTION########################################
##################################################################################################


########################### Preparing scan file ###############################

	#Listing all files in the directory (Networks,genes_results and mutload results)
	if (is.null(arg$network)) {
	networks<-gsub('.{5}$', '', list.files(pattern=paste0(".json")))
	
	} else { networks<-arg$network	}
	
	
	#This commented out part can be used to exclude analysis of file results that exist in folder
	#genes_results_files<-list.files(pattern=paste0(".*_genes_results"))
	#mutload_results_files<-list.files(pattern=paste0(".*_mutload_results"))

	# Inferring network files that have not been processed yet
	#networks_to_process_genes<-sapply(networks,function (x) { sum(grepl(pattern = x,genes_results_files))} )
	#networks_to_process_genes<-names(networks_to_process_genes[networks_to_process_genes==0])
	#networks_to_process_mutload<-sapply(networks,function (x) { sum(grepl(pattern = x,mutload_results_files))} )
	#networks_to_process_mutload<-names(networks_to_process_mutload[networks_to_process_mutload==0])
	#scan$genes_connectivity<-scan$networks %in% networks_to_process_genes
	#scan$mutload_connectivity<-scan$networks %in% networks_to_process_mutload
	#scan<-read.csv("dict.csv",as.is=T)
	
	scan<-data.frame(networks=networks,resolution=NA,gain=NA,original_samples=NA,first_connected_samples=NA,samples_threshold=NA,p_0.05=NA,q_0.1=NA,q_0.15=NA,q_0.2=NA,mutload=NA,stringsAsFactors = F)

	scan$resolution<-as.numeric(sapply(scan$networks, function (x) {
	  strsplit(x,"_")[[1]][4]
	}))

	scan$gain<-as.numeric(sapply(scan$networks, function (x) {
	  strsplit(x,"_")[[1]][5]
	}))

	if (arg$scan!=0) { #Removing network files file for test mode
	  scan<-scan[1:arg$scan,]
	}









############################################################################################################
################### Running connectivity analysis for networks in scan table ##############################
############################################################################################################
	
count<-0

#for (file in scan$networks[scan$mutload_connectivity]) {
for (file in scan$networks) {
	   count<-count+1  
	   print("*********************************************")
	   print (paste("Analyzing network:",file,"-",count,"out of",nrow(scan)))
	   print (paste("Number of permutations:",arg$permutations))
	   print (paste("Number of CPU cores:",arg$cores))
	  
	  #run_line<-paste("Rscript", arg$connectivity, "-p",arg$permutations,"-h TRUE -n",file,"-m",arg$matrix,"-q",arg$cores,"-k",arg$chunk)
	  #system(run_line)

	 #Parsing and loading, gexf(edge file) and json (nodes file) to memory.
	 print ("Parsing json and gexf files")
	 graph_name<-file
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
	 matrix1<-mat_non_syn_bin[samples_relabling_table[,1],] #Subseting matrix to contain only samples in first connected graph 

	 #Extracting columns from arguments
	 columns<-column_range(arg$columns)
	 #Removing columns below samples_threshold from the first connected graph
	 matrix1<-matrix1[,columns,drop=FALSE] #Subsetting for selected columns



	 #Choosing genes based on score
	 samples_of_interest<-rownames(matrix1)

	 # Taking record of sample sizes
	
	 scan[scan$networks==file,]$original_samples<-length(all_samples)
	 scan[scan$networks==file,]$first_connected_samples<-length(samples_of_interest)

	 ###############################################
	 ##################RESCALING WAS HERE####################
	 ###############################################

	
	
	
	##########################################
	############THRESHOLDING SECTION###########
    ##########################################


	# selecting genes based on thresholds
	
     samples_threshold<-ceiling(arg$samples_threshold*length(samples_of_interest))
     scan[scan$networks==file,]$samples_threshold<-samples_threshold
     print (paste("Samples in original dataset:",length(all_samples)))
     print (paste("Samples in first connected graph:",length(samples_of_interest)))
     print (paste("Samples threshold is set to:",samples_threshold))
     print (paste("Top genes score threshold is set to:",arg$g_score_threshold))
     
	 genes_number_of_samples<-apply(matrix1,2,function (x) sum(x!=0)) #Counting non_zero samples for each column
	 genes_below_samples_threshold<-names(which(genes_number_of_samples<samples_threshold))
	 genes_above_samples_threshold<-names(which(genes_number_of_samples>=samples_threshold)) #For filtering by number  of mutations exist in a sample


	 

	 if (arg$score_type=="lam") {
	 g_score<-g_score_calc(arg$score_type,samples_of_interest,all_genes) #all_genes_Lambda scores needs all genes into account 
	 } else
	 g_score<-g_score_calc(arg$score_type,samples_of_interest,genes_above_samples_threshold) #Syn/old only over sample thresholded genes 

	 columns_of_interest<-head(sort(g_score,decreasing = T),arg$g_score_threshold) #Filtering by g-score

	 print(paste0("Columns above threshold: ",length(columns_of_interest)))

	 
	 
	 
	 matrix1<-matrix1[,names(columns_of_interest),drop=FALSE] #Subsetting matrix to have above threshold columns
     
     if (arg$mutload==TRUE) {
     #Adding to matrix1 a column with mutation rate, this will be used to assess mutload mutated samples. 
	 	
		  matrix1<-as.matrix(mutload_dist[samples_of_interest],drop=FALSE)
		  colnames(matrix1)<-"mutLoad"
		  columns_of_interest<-"mutLoad"
	 }
		
	 #Initializing results file name and unique id
	 unique_id<-round(runif(1, min = 111111, max = 222222),0)
	 file_prefix<-paste0(file,"_",arg$matrix,"-",unique_id,"-",Sys.Date())
	 print(paste("File unique identifier:",unique_id))

	 

	 if (arg$mutload==12) {
	 print ("CCCCC")
	 #Adding to matrix1 a column with mutation rate, this will be used to assess mutload mutated samples. 
     
		 
		  mutLoad<-mat_non_syn+mat_syn #Total number of point mutations
			  if (arg$rescale!=0) {
		   png(paste0(file_prefix,"_mutLoad_Rescaled.png"))
		   invisible(dev.off())  
		  } else {
		   png(paste0(file_prefix,"_mutLoad_NoRescaling.png"))
		   hist(log10(mutLoad),breaks = 100,main="Before rescaling")
		   invisible(dev.off())
		  }
		  
		  mutLoad<-rowSums(mutLoad)[samples_of_interest] #rownames matrix1 is important to account only for samples_of_interes
		  matrix1<-as.matrix(mutLoad,drop=FALSE)
		  colnames(matrix1)<-"mutLoad"
		  columns_of_interest<-"mutLoad"
		  
	 }

	# info_cols was here

	 #Printing thresholded genes
	 #thresholded_genes1<-genes_below_samples_threshold
	 #thresholded_genes2<-setdiff(genes_above_samples_threshold,names(columns_of_interest))
	 #write.csv(thresholded_genes1,paste0(file_prefix,"_thresholded_genes_samples.csv"))
	 #write.csv(thresholded_genes2,paste0(file_prefix,"_thresholded_genes_score.csv"))



	 
	 #Writing log file
	 suppressWarnings(write.table(as.character(arg) ,paste0(file_prefix,"_log.csv"),append=TRUE))
	 suppressWarnings(write.table(paste("Number of permutations: ",arg$permutations),paste0(file_prefix,"_log.csv"),append=TRUE))
	 suppressWarnings(write.table(paste("Samples threshold: ",arg$samples_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
	 suppressWarnings(write.table(paste("g_score threshold: ",arg$g_score_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
	 suppressWarnings(write.table(paste("Columns above threshold:",length(columns_of_interest)),paste0(file_prefix,"_log.csv"),append=TRUE))
	 suppressWarnings(write.table(paste0("Original sample size:",length(all_samples)),paste0(file_prefix,"_log.csv"),append=TRUE))
	 suppressWarnings(write.table(paste0("First connected sample size:",length(samples_of_interest)),paste0(file_prefix,"_log.csv"),append=TRUE))




	 #logger <- create.logger(logfile = 'debugging.log', level = 1)
	 #info(logger,paste("Number of permutations: ",arg$permutations))
	 #info(logger,paste("Samples threshold: ",arg$samples_threshold))


	 permutations<-arg$permutations
	 edges1<-edges[,1] #Nodes i
	 edges2<-edges[,2] #Nodes j
	 num_nodes<-length(nodes)


	 print (Sys.time())
	 ptm<-proc.time()




	 ######################### CONNECTIVITY ANALYSIS ##############################################3
	 
	 print ("Starting connectivity analysis:")
	 if (arg$mutload==TRUE) {
		print (paste("Mutload Connectivity for Graph"))
	 	ans<-connectivity_analysis(columns_of_interest,matrix1)
	 	final_results<-results_file(ans)
	 	scan[scan$networks==file,]$mutload<-final_results[,"p_value"]

	 
	 } else {   # Genes analysys
			print (paste("Genes Connectivity for Graph"))
			ans<-connectivity_analysis(columns_of_interest,matrix1)
	 		final_results<-results_file(ans)
	 		scan[scan$networks==file,]$p_0.05<-sum(final_results[,"p_value"]<=0.05)
	 		scan[scan$networks==file,]$q_0.1<-sum(final_results[,"q_value"]<=0.1)
	 		scan[scan$networks==file,]$q_0.15<-sum(final_results[,"q_value"]<=0.15)
	 		scan[scan$networks==file,]$q_0.2<-sum(final_results[,"q_value"]<=0.2)
	 }
	 			
	 

		################ pii_value generations goes here - pushed to end##################
	



	 
	 ##################Synonymous control###############
	 if (arg$syn_control==TRUE & length(columns_of_interest)!=0) {
	 print ("Starting control connectivity analysis:")
	 matrix1<-ifelse(mat_syn>0,1,0) #Using synonymous matrix as reference
	 matrix1<-matrix1[samples_of_interest,names(columns_of_interest),drop=FALSE] # Subsetting for samples of interest
	 ans<-connectivity_analysis(columns_of_interest,matrix1) #Running connectivity analysis
	 final_results_control<-results_file(ans)

	 #Coercing non_syn and control results
	 final_results_control<-final_results_control[,c("n_samples","p_value","q_value"),drop=FALSE]
	 colnames(final_results_control)<-c("n_samples_con","p_value_con","q_value_con")
	 missing_genes<-setdiff(rownames(final_results),rownames(final_results_control)) #Genes that do not exist in final_Results needed to be completed with NA and zeros

	 missing_n_samples_con<-colSums(matrix1[,missing_genes])
	 missing_p_value_con<-rep(NA,length(missing_genes))
	 missing_q_value_con<-rep(NA,length(missing_genes))
	 missing_x<-data.frame(missing_n_samples_con,missing_p_value_con,missing_q_value_con)
	 final_results_control<-rbind(final_results_control,as.matrix(missing_x))
	 final_results_control<-final_results_control[rownames(final_results),,drop=FALSE]
	 final_results<-cbind(final_results,final_results_control)
	 final_results<-final_results[,c("Gene_Symbol","EntrezID","c_value","p_value","n_samples","q_value","p_value_con","n_samples_con","g_score_syn"),drop=FALSE]

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
	 speed_index<-run_t[3]*500/arg$permutations/length(columns)
	 print(paste("Runtime in seconds:",run_t[3]))
	 print(paste("Speed index (calc time for 500 permutations):",speed_index)) 
	    
}



#Writing scanner summary file:
write.csv(scan,"scan_summary.csv")


######################################## Plotting section######################################



#Number of samples per graph plot:
ggplot(scan, aes(x=resolution, y=gain, label=first_connected_samples)) + 
  #scale_color_gradient2(low = 'white', mid='yellow', high = 'red') +
  geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle(paste("Original number of samples:",scan$original_samples[1])) +
  ggsave(filename = paste0("First_component_samples.png"))  






if (arg$mutload==FALSE) {  #Connectivity plots and number_of_Events


	 #connectivity_q_value_plot
	 q_threshold_range<-c(0.1,0.15,0.2)
	 for (threshold in q_threshold_range) {

		q_value_dist<-scan[,paste0("q_",threshold)]
		title<-paste("Genes_results_q_value <=",threshold, "Permutations=",arg$permutations)

		ggplot(scan, aes(x=resolution, y=gain, color=q_value_dist, label=q_value_dist)) + 
		scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
		geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle(title) +
		ggsave(filename = paste0("Genes_results_q_value","_",threshold,".png"))   
 
	 }

	 #connectivity_p_value_plot
	 p_value_dist<-scan$p_0.05
	 ggplot(scan, aes(x=resolution, y=gain, color=p_value_dist, label=p_value_dist)) + 
	  scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
	  geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle("Genes_results_p_value<=0.05") +
	  ggsave(filename = paste0("Genes_results_p_value_0.05.png"))    




	 genes_results_files<-
	  sapply(scan$networks,function (x) {
		results<-list.files(pattern=paste0("^",x,".*_genes_results"))
		if (length(results)==0) {results<-NA}
		return(results)
	  })




	 ################ Number of events per gene###########################

	 #Calculates sum of a particular feature for a list of genes across scan results
	 number_of_events<-function(genes_results_files,feature,threshold) {

	  genes<-sapply(genes_results_files,function (file)
	  {
		if (!is.na(file)) {
		  results<-read.csv(file,as.is=T)
		  results<-results$Gene_Symbol[results[,feature]<=threshold]  
		}

	  })

	  genes<-unlist(genes)
	  genes<-genes[complete.cases(genes)]
	  unique_genes<-unique(genes)

	  genes_events<-sapply(unique_genes,function (x) sum(x==genes))
	  return(sort(genes_events,decreasing = T))

	 }

	 #Generating number of events summary file
	 q_value_0.1<-number_of_events(genes_results_files,"q_value",0.1)
	 q_value_0.15<-number_of_events(genes_results_files,"q_value",0.15)
	 q_value_0.2<-number_of_events(genes_results_files,"q_value",0.2)
	 p_value_0.05<-number_of_events(genes_results_files,"p_value",0.05)

	 n<-max(length(q_value_0.1),length(q_value_0.15),length(q_value_0.2),length(p_value_0.05))
	 length(q_value_0.1)<-n ; length(q_value_0.15) <-n; length(q_value_0.2) <-n; length(p_value_0.05) <-n

	 events<-data.frame(
			  q_value_0.1,
			  q_value_0.15,
			  q_value_0.2,
			  p_value_0.05
			  )

	 colnames(events)<-c("q_value_0.1","q_value_0.15","q_value_0.2","p_value_0.05")
	 write.csv(events,"number_of_events.csv")



	    
} else {   # MUTLOAD PLOT
	ggplot(scan, aes(x=factor(resolution), y=gain)) + 
  	geom_point(size=5,aes(color=scan$mutload<=0.05)) + geom_text(label=scan$mutload,vjust=1.6)+theme_bw() + ggtitle("Mutational Load Connectivity") + 
  	guides(color = guide_legend(title = paste("mutload <= 0.05"),
                title.theme = element_text(size=10,angle=0,color="blue"))) +  ggsave(filename = "mutload_grid.png")

}




###############MOVING FILES TO Results DIR####################



print ("Moving files to Results Directory")
if (!file.exists("Results")) {dir.create("Results")}
if (file.exists("Rplots.pdf")) {file.remove("Rplots.pdf")}

csv_files<-list.files(pattern = "*.csv")
png_files<-list.files(pattern = "*.png")

files_to_move_to_results<-c(csv_files,png_files)
x<-file.rename(files_to_move_to_results,paste0("Results/",files_to_move_to_results))
if (sum(x)==length(files_to_move_to_results)) {
  print ("All results files moved to Results dir, archiving files")
  tar(paste0("Results_",arg$matrix,".tar.gz"),"Results")
} else {
        print ("This files were not moved to Results dir:")
        print (files_to_move_to_results[!x])
}




#####################OLD FUNCTIONS AND PROCEDURES########################
##########################################################################

 #############################################################
	 #Generating pii_values table
	 #pi_values_table<-NULL
	 #for (i in 1:length(ans))
	 #  pi_values_table<-cbind(pi_values_table,ans[[i]][[2]]) #Extracting pii_values from ans
	 #colnames(pi_values_table)<-pi_values_table[1,]
	 #pi_values_table<-pi_values_table[,!pi_zero_genes,drop=FALSE]
	 #############################################################

	 #Writing final results and pii_values files
	 #write.table(final_results,paste0(file_prefix,"_results_final.csv"),row.names=FALSE,sep=",")
	 #write.table(pi_values_table,paste0(file_prefix,"_pii_values.csv"),sep=",",row.names=FALSE,col.names=FALSE)







#extract_value<-function (files,feature,threshold)  {
  #Gets a feature(p_value) and a threshold (0.05) and extract the number of observations across list of files below that threshold for this feature
#  ans<-sapply(files,function (file) {
#    if (length(file)!=1 | is.na(file)) {
#      results<-NA
#    } else {
#      results<-read.csv(file,as.is=T)[,feature]
#      results<-sum(results<=threshold,na.rm=T)
#    }
#  })
#  return (as.numeric(ans))
#}



#genes q_value plot
#threshold_range<-c(0.1,0.15,0.2)
#for (threshold in threshold_range) {
  
#  q_value_dist<-extract_value(genes_results_files,"q_value",threshold)
#  title<-paste("Genes_results_q_value <=",threshold, "Permutations=",arg$permutations)
#  ggplot(scan, aes(x=resolution, y=gain, color=q_value_dist, label=q_value_dist)) + 
#    scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
#    geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle(title) +
#    ggsave(filename = paste0("Genes_results_q_value","_",threshold,".png"))    
#}







