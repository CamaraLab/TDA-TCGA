

column_range<-function(col_range,matrix)
  #Gets column range from arg$column and parse it. Also removes columns below threshold
{
  
  if (col_range=="all") #Check if all is supplied
  {
    x<-1:ncol(matrix)
    
  } else if (grepl(":",col_range)==TRUE) { #If range x:x is supplied
    
    x<-as.numeric(strsplit(col_range,":")[[1]][1]:strsplit(col_range,":")[[1]][2])
    
  } else if (suppressWarnings(!is.na(as.numeric(col_range)))) { # If column is all numbers  
    x<-as.numeric(col_range)  
  } else if (!is.na(match(col_range,colnames(matrix)))) { #If Gene name exist, convert to column number
    x<-which(colnames(matrix)==col_range)       
  } else stop ("Column name not found")
  
  
  
  #Validating column in range:
  if (max(x)>ncol(matrix)) 
  {
    stop (paste0("Columns out of range"," Available range: 1:",ncol(matrix)))
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
  #c<-sum(pi_column[edges1]*pi_column[edges2]))
  
  c_vector<-c_vector*(num_nodes/(num_nodes-1))
} 




connectivity_analysis<-function(columns_of_interest,matrix) {
  
  print (paste("Preparing parallel environment. Acquiring",arg$cores,"Cores"))
  cl <- makeCluster(as.numeric(arg$cores))
  varlist=c("file_prefix","c_calc_fast","c_calc_fast","pii_matrix","e_matrix","edges1","edges2","samples","permutations","num_nodes","perm_values","arg","nodes","matrix","info_cols","columns")
  clusterExport(cl=cl, varlist=varlist,envir=environment())
  
  columns<-seq_along(columns_of_interest) #Subsetting columns
  split.column<-split(columns,ceiling(seq_along(columns)/arg$chunk))
  
  #count<-0
  #ans<-parLapply(cl,split.column,function (columns_range)  {
  ans<-lapply(split.column,function (columns_range)  {
    #calculating c-scores and p-values for each chunk of columns
    #count<-count+1
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
    
    
    if (1==1) {
      e_list<-perm_values_list
    } else e_list<-lapply(perm_values_list,function(x) e_matrix(nodes,as.matrix(x))) #columns_range elements in the list. Each element is a matrix representing pi_values of a gene. rows are nodes, columns are permutations. first column is non permuted. 
    
    pi_list<-lapply(e_list,function(x) pii_matrix(as.matrix(x))) #columns_range elements in the list. Each element is a matrix representing pi_values of a gene. rows are nodes, columns are permutations. first column is non permuted.
    c_vec_list<-lapply(pi_list,function (pi_matrix) c_calc_fast(as.matrix(pi_matrix)))
    
    e_mean<-sapply(e_list,function (x) mean(x[,1])) #Taking mean of the first column (not permutations)
    e_sd<-sapply(e_list,function (x) sd(x[,1]))
    pi_values<-sapply(pi_list,function (x) x[,1]) #Is a matrix,each row is a node, each column in the matrix is pi values of a gene across nodes.
    pi_frac<-apply(pi_values,2,function (x) sum(x!=0)/length(x))
    n_samples<-apply(matrix2,2,function (x) sum(x!=0))
    c_value<-sapply(c_vec_list,function (c_vec) c_vec[1]) #the first position is the connectivity value, the later are c_value for each permutation
    p_value<-sapply(c_vec_list,function(c_vec) {
      p_value<-sum(c_vec>c_vec[1])/permutations})
    Genes<-colnames(matrix2)  
    
    output<-cbind(Genes,c_value,p_value,pi_frac,n_samples,e_mean,e_sd) #The variable names should match info_cols
    output1<-t(as.data.frame(c_vec_list))
    rownames(output1)<-Genes
    
    write.csv(output1,paste0("output_",columns_range[1],".csv"))
    #write.csv(output1,"output1.csv")
    #write.csv("output1.csv",output1)
    #write.table(output,paste0(file_prefix,"_results_rolling.csv"),append=TRUE,sep=",",col.names=FALSE,row.names=FALSE)
    
    #Ans is columns_range length list. Each element contain to variables.
    #First variabl is "output" which is results matrix. the second element is pii_values matrix, its columns are genes and rows are nodes 
    
    ans<-list(output,rbind(Genes,pi_values),output1)
    
    return(ans) 
    
  })
  print("Releasing cores")
  stopCluster(cl) # Releasing acquired CPU cores
  return(ans)
}

#x<-ans[[1]][[3]]
#for (i in 2:length(ans)) {
##  print(i)
  #dim(ans[[i]][[3]])
  #x<-rbind(x,ans[[i]][[3]])
#}
#write.csv(x,"x.csv")



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



epsilon_p_value<-function(results_files) {
  #Takes in a list of results files and returns average p_values of all the genes
  Gene_Symbol<-sapply(strsplit(names(columns_of_interest),"|",fixed = TRUE),"[[",1)
  Gene_Symbol<-substring(Gene_Symbol,5)
  p_list<-lapply(Gene_Symbol,function(x) NULL)
  names(p_list)<-Gene_Symbol
  
  for (i in results_files) {
    results<-read.csv(i)
    for (gene in names(p_list)) {
      p_list[[gene]]<-c(p_list[[gene]],results$p_value[results$Gene_Symbol==gene])
    } 
  }
  
  ans<-p_list
  
}


epsilon_q_value<-function(results_files) {
  #Takes in a list of results files and returns average p_values of all the genes
  Gene_Symbol<-sapply(strsplit(names(columns_of_interest),"|",fixed = TRUE),"[[",1)
  Gene_Symbol<-substring(Gene_Symbol,5)
  p_list<-lapply(Gene_Symbol,function(x) NULL)
  names(p_list)<-Gene_Symbol
  
  for (i in results_files) {
    results<-read.csv(i)
    for (gene in names(p_list)) {
      p_list[[gene]]<-c(p_list[[gene]],results$q_value[results$Gene_Symbol==gene])
    } 
  }
  
  ans<-p_list
  
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
setwd("c:/Users/Udi/SkyDrive/TCGA_CURATED/Rips/")
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
arg<-list(20,NULL,NULL,NULL,NULL,"STADTRIM.h5","all",10,detectCores(),FALSE,TRUE,NULL,0.06,100,"syn","Annotations.csv",FALSE,FALSE,0,"PROCESSED_COAD_hgsc.bcm.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2_OCT_16_2015.maf")
names(arg)<-c("epsilon","cut","topgenes","network","scan","matrix","columns","permutations","cores","log2","fdr","chunk","samples_threshold","g_score_threshold","score_type","anno","mutload","syn_control","rescale","maf")

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


#amples_to_remove1<-c("TCGA-FP-8209-01","TCGA-FP-8210-01","TCGA-BR-A4IZ-01")
#samples_to_remove2<-c("TCGA-BR-6452-01","TCGA-BR-8680-01","TCGA-CG-5721-01","TCGA-BR-8487-01","TCGA-BR-4361-01")
#samples_to_remove<-c("TCGA-BR-6452-01","TCGA-BR-8680-01","TCGA-CG-5721-01","TCGA-BR-8487-01","TCGA-BR-4361-01")
#r<-match(samples_to_remove,all_samples)
#mat_tpm<-mat_tpm[-r,]
#mat_syn<-mat_syn[-r,]
#mat_non_syn<-mat_non_syn[-r,]
#mat_non_syn_bin<-mat_non_syn_bin[-r,]

#all_samples<-all_samples[-r]

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
info_cols<-t(c("Genes","c_value","p_value","pi_frac","n_samples","e_mean","e_sd")) 






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

exp_var<-apply(mat_tpm[,exp_mean>arg$cut],2,som)
top5000<-names(head(sort(round(exp_var,3),decreasing = T),arg$topgenes))

exp_top5000<-mat_tpm[,top5000]
exp_top5000<-t(exp_top5000)
#cor_exp_top5000<-1-cor(exp_top5000)
if (arg$dist==1) {
	cor_exp_top5000<-1-cor(exp_top5000,method="pearson")
} 

if (arg$dist==2) {
	cor_exp_top5000<-1-cor(exp_top5000,method="spearman")
} 

if (arg$dist==3) {
	cor_exp_top5000<-as.matrix(dist(t(exp_top5000),upper = T))
}

#require("hopach")
#cor_exp_top5000<-1-as.matrix(mutualInfo(t(exp_top5000)))
  #as.matrix(distancematrix(t(exp_top5000),"eucli"))
#g<-1-cor(exp_top5000)
#sum(is.na(cor_exp_top5000))

epsilon_min<-min(cor_exp_top5000)
epsilon_max<-max(cor_exp_top5000)
epsilon_set<-round(seq(epsilon_min,epsilon_max,length.out=arg$epsilon),2)
print (epsilon_set)
scan<-data.frame(networks=epsilon_set,resolution=NA,gain=NA,original_samples=NA,first_connected_samples=NA,edges_num=NA,samples_threshold=NA,above_samples_threshold=NA,above_gscore_threshold=NA,p_0.05=NA,q_0.1=NA,q_0.15=NA,q_0.2=NA,mutload=NA,parsing_time=NA,connectivity_time=NA,uid=NA,stringsAsFactors = F)
  
if (arg$scan!=0) { #Removing network files file for test mode
  scan<-scan[1:arg$scan,]
}





############################################################################################################
################### Running connectivity analysis for networks in scan table ##############################
############################################################################################################

count<-0
c_matrix_list<-list()
#for (file in scan$networks[scan$mutload_connectivity]) {
for (file in scan$networks) {
  count<-count+1  
  print("*********************************************")
  print (Sys.time())
  print (paste("Analyzing network:",file,"-",count,"out of",nrow(scan)))
  print (paste("Number of permutations:",arg$permutations))
  print (paste("Number of CPU cores:",arg$cores))
  
  adj_mat<-ifelse(cor_exp_top5000<=file,1,0)
  adj_mat[lower.tri(adj_mat,diag = TRUE)]<-0
  rownames(adj_mat)<-1:nrow(adj_mat)
  colnames(adj_mat)<-1:ncol(adj_mat)
  graph_igraph<-graph.adjacency(adj_mat,diag = F)
  edges<-get.edgelist(graph_igraph,names=FALSE)
  nodes<-lapply(1:nrow(adj_mat),function (x) x)
  names(nodes)<-1:nrow(adj_mat)
  
  edges_num<-nrow(edges)
  print (paste("Number of edges in network:",edges_num))
  scan[scan$networks==file,]$edges_num<-edges_num
  
  samples<-unique(unlist(nodes))
  matrix1<-mat_non_syn_bin[samples,,drop=FALSE] #Subseting matrix to contain only samples in first connected graph 
  
  #Extracting columns from arguments
  columns<-column_range(arg$columns,matrix1)
  #Removing columns below samples_threshold from the first connected graph
  matrix1<-matrix1[,columns,drop=FALSE] #Subsetting for selected columns
    #Choosing genes based on score
  samples_of_interest<-rownames(matrix1)
  
  # Taking record of sample sizes
  
  scan[scan$networks==file,]$original_samples<-length(all_samples)
  scan[scan$networks==file,]$first_connected_samples<-length(samples_of_interest)
  
  
  
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
  
  
    
  g_score<-g_score_calc(arg$score_type,samples_of_interest,genes_above_samples_threshold) #Syn/old only over sample thresholded genes
  columns_of_interest<-head(sort(g_score,decreasing = T),arg$g_score_threshold) #Filtering by g-score
  
  print(paste0("Columns above threshold: ",length(columns_of_interest)))
  
  scan[scan$networks==file,]$above_samples_threshold<-length(genes_above_samples_threshold)
  scan[scan$networks==file,]$above_gscore_threshold<-length(columns_of_interest)
  
  
  matrix1<-matrix1[,names(columns_of_interest),drop=FALSE] #Subsetting matrix to have above threshold columns
  
  if (arg$mutload==TRUE) {
    #Adding to matrix1 a column with mutation rate, this will be used to assess mutload mutated samples. 
    
    matrix1<-as.matrix(mutload_dist[samples_of_interest],drop=FALSE)
    colnames(matrix1)<-"mutLoad"
    columns_of_interest<-"mutLoad"
  }
  
  #Initializing results file name and unique id
  unique_id<-round(runif(1, min = 111111, max = 222222),0)
  file_prefix<-paste0(file,"_",PROJECT_NAME,"-",unique_id,"_",guid,"-",Sys.Date())
  print(paste("File unique identifier:",unique_id))
  scan[scan$networks==file,]$uid<-unique_id
  
  
  
    
  #Writing log file
  suppressWarnings(write.table(as.character(arg) ,paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("Number of permutations: ",arg$permutations),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("Samples threshold: ",arg$samples_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("g_score threshold: ",arg$g_score_threshold),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste("Columns above threshold:",length(columns_of_interest)),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste0("Original sample size:",length(all_samples)),paste0(file_prefix,"_log.csv"),append=TRUE))
  suppressWarnings(write.table(paste0("First connected sample size:",length(samples_of_interest)),paste0(file_prefix,"_log.csv"),append=TRUE))
  
  
  
  permutations<-arg$permutations
  edges1<-edges[,1] #Nodes i
  edges2<-edges[,2] #Nodes j
  num_nodes<-length(nodes)
  
  
  
  ptm<-proc.time()
  
  
  
  
  ######################### CONNECTIVITY ANALYSIS ##############################################3
  
  if (edges_num!=0) {
    
    
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
    
    
    c_matrix<-ans[[1]][[3]]
    for (i in 2:length(ans)) {
      ##  print(i)
      #dim(ans[[i]][[3]])
      c_matrix<-rbind(c_matrix,ans[[i]][[3]])
    }
    c_matrix_file<-paste0(file_prefix,"_c_matrix.csv")
    write.csv(c_matrix,c_matrix_file)
    
    ################ pii_value generations goes here - pushed to end##################
    c_matrix_list[[as.character(file)]]<-c_matrix    
    
    
    
    ##################Synonymous control###############
    if (arg$syn_control==TRUE & length(columns_of_interest)!=0) {
      print ("Starting control connectivity analysis:")
      matrix1<-ifelse(mat_syn>0,1,0) #Using synonymous matrix as reference
      matrix1<-matrix1[samples_of_interest,names(columns_of_interest),drop=FALSE] # Subsetting for samples of interest
      ans<-connectivity_analysis(columns_of_interest,matrix1) #Running connectivity analysis
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
    scan[scan$networks==file,]$connectivity_time<-run_t[3]
    speed_index<-run_t[3]*500/arg$permutations/length(columns)
    print(paste("Runtime in seconds:",run_t[3]))
    print(paste("Speed index (calc time for 500 permutations):",speed_index)) 
    
  } else {
    print ("NO EDGES IN THE GRAPH ONLY SCAN FILE IS PRODUCED")
    scan[scan$networks==file,][,7:ncol(scan)]<-0  # NO EDGES IN THE GRAPH
    
  }
  
  
  
}



#Writing scanner summary file:
write.csv(scan,paste0("scan_summary_",guid,".csv"))


######################################## Plotting section######################################


if (arg$mutload==FALSE) {  #Connectivity plots and number_of_Events
  
    
genes_results_files<-list.files(pattern=paste0(scan$uid,".*_genes_results.csv"))


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
  
  gene_order<-names(q_value_0.1)
  
  networks_with_no_edges<-scan$networks[scan$edges_num==0] 
  ee<-setdiff(epsilon_set,networks_with_no_edges)


  dp<-epsilon_p_value(genes_results_files)[gene_order]
  dq<-epsilon_q_value(genes_results_files)[gene_order]
  average_p<-sapply(dp,mean,na.rm=T)
  average_q<-sapply(dq,mean,na.rm=T)
  q_average_p<-p.adjust(average_p,method = "fdr")
  
  epsilon_dist_p<-as.data.frame(t(as.data.frame(dp)))
  epsilon_dist_q<-as.data.frame(t(as.data.frame(dq)))
  colnames(epsilon_dist_p)<-ee
  colnames(epsilon_dist_q)<-ee
  epsilon_dist_p$average_p<-average_p
  epsilon_dist_q$average_q<-average_q


  write.csv(epsilon_dist_p,paste0("epsilon_dist_p_",guid,".csv"))
  write.csv(epsilon_dist_q,paste0("epsilon_dist_q_",guid,".csv"))

  
  k_som<-epsilon_dist_p
  k_som<-as.matrix(k_som)
  c<-melt(k_som)
  colnames(c)<-c("gene","epsilon","p_value")

  c$plot<-c$p_value<=0.05
  c<-c[c$plot,]

  
  
x<-ncol(epsilon_dist_p)
cut_left<-1+ceiling(x/10)
cut_right<-(x-(cut_left-1)*2)
  
  ggplot(c,aes(x=epsilon,y=gene)) + geom_point() + 
  geom_vline(xintercept = c((cut_left-0.2),(cut_right+0.2)),col="red") +
  ggsave(filename=paste0("epsilon_plot_p_,",guid,".png"))


k_som<-epsilon_dist_q
k_som<-as.matrix(k_som)
c<-melt(k_som)
colnames(c)<-c("gene","epsilon","q_value")

c$plot<-c$q_value<=0.15
c<-c[c$plot,]

x<-ncol(epsilon_dist_q)
cut_left<-1+ceiling(x/10)
cut_right<-(x-(cut_left-1)*2)


ggplot(c,aes(x=epsilon,y=gene)) + geom_point() +
  geom_vline(xintercept = c((cut_left),(cut_right)),col="red") +
  ggsave(filename = paste0("epsilon_plot_q_",guid,".png"))

#PLot after removing outlier epsilons

k_som<-epsilon_dist_q[cut_left:cut_right]
k_som<-as.matrix(k_som)
c<-melt(k_som)
colnames(c)<-c("gene","epsilon","q_value")


write.csv(k_som,paste0("epsilon_dist_q_cut_",guid,".csv"))
write.csv(epsilon_dist_q,paste0("epsilon_dist_q_",guid,".csv"))

k_som<-as.data.frame(k_som)
#k_som<-as.matrix(k_som)
k_som$average_q<-rowMeans(k_som)
k_som$frequency<-apply(k_som[,1:(ncol(k_som)-1)],1,function (x) sum(x<=0.15))
write.csv(k_som,paste0("epsilon_dist_q_cut_",guid,".csv"))



c$plot<-c$q_value<=0.15
c<-c[c$plot,]
ggplot(c,aes(x=epsilon,y=gene)) + geom_point() +ggsave(filename = paste0("epsilon_plot_q_cut_",guid,".png"))


#plot p value after cut

k_som<-epsilon_dist_p[cut_left:cut_right]
k_som<-as.matrix(k_som)
c<-melt(k_som)
colnames(c)<-c("gene","epsilon","p_value")


k_som<-as.data.frame(k_som)
k_som$average_p<-rowMeans(k_som)
k_som$frequency<-apply(k_som[,1:(ncol(k_som)-1)],1,function (x) sum(x<=0.05))
write.csv(k_som,paste0("epsilon_dist_p_cut_",guid,".csv"))



c$plot<-c$p_value<=0.05
c<-c[c$plot,]
ggplot(c,aes(x=epsilon,y=gene)) + geom_point() +ggsave(filename = paste0("epsilon_plot_p_cut_",guid,".png"))







  n<-max(length(q_value_0.1),length(q_value_0.15),length(q_value_0.2),length(p_value_0.05))
  length(q_value_0.1)<-n ; length(q_value_0.15) <-n; length(q_value_0.2) <-n; length(p_value_0.05) <-n
  
  events<-data.frame(
    q_value_0.1,
    q_value_0.15,
    q_value_0.2,
    p_value_0.05
  )
  
  colnames(events)<-c("q_value_0.1","q_value_0.15","q_value_0.2","p_value_0.05")
  write.csv(events,paste0("number_of_events_",guid,".csv"))
  

  
  
  
} else {   # MUTLOAD PLOT
  print(as.numeric(scan$mutload))
  scan$mutload<-scan$mutload
  ggplot(scan, aes(x=networks, y=mutload)) + 
    geom_point(size=5,aes(color=scan$mutload<=0.05)) + geom_text(label=scan$mutload,vjust=1.6)+theme_bw() + ggtitle("Mutational Load Connectivity") + 
    guides(color = guide_legend(title = paste("mutload <= 0.05"),
                                title.theme = element_text(size=10,angle=0,color="blue"))) +  ggsave(filename = paste0("mutload_grid_",guid,".png"))
  
}




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

