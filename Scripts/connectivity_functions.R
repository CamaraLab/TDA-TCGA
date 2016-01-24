

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
    c_value<-sapply(c_vec_list,function (c_vec) c_vec[1])
    p_value<-sapply(c_vec_list,function(c_vec) {
      p_value<-sum(c_vec>c_vec[1])/permutations})
    Genes<-colnames(matrix2)  
    
    output<-cbind(Genes,c_value,p_value,pi_frac,n_samples,e_mean,e_sd) #The variable names should match info_cols
    #write.csv("output1.csv",outertput1)
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



