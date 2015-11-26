#rm(list=ls())
#setwd("c:/Users/Udi/SkyDrive/TCGA_CURATED/Rips/")
#del<-list.files(pattern = "*.csv")
#unlink(del)
source('C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Scripts/connectivity_functions.R')
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
arg<-list(10,NULL,NULL,NULL,NULL,"STADCUR.h5","all",5,detectCores(),FALSE,TRUE,NULL,0.06,100,"syn","Annotations.csv",FALSE,FALSE,2,"PROCESSED_COAD_hgsc.bcm.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2_OCT_16_2015.maf")
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


samples_to_remove1<-c("TCGA-FP-8209-01","TCGA-FP-8210-01","TCGA-BR-A4IZ-01")
samples_to_remove2<-c("TCGA-BR-6452-01","TCGA-BR-8680-01","TCGA-CG-5721-01","TCGA-BR-8487-01","TCGA-BR-4361-01")
samples_to_remove<-c("TCGA-BR-6452-01","TCGA-BR-8680-01","TCGA-CG-5721-01","TCGA-BR-8487-01","TCGA-BR-4361-01")
r<-match(samples_to_remove,all_samples)
mat_tpm<-mat_tpm[-r,]
mat_syn<-mat_syn[-r,]
mat_non_syn<-mat_non_syn[-r,]
mat_non_syn_bin<-mat_non_syn_bin[-r,]

all_samples<-all_samples[-r]

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
	cor_exp_top5000<-as.matrix(dist(exp_top5000,upper = T))
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
    
    
    
    ################ pii_value generations goes here - pushed to end##################
    
    
    
    
    
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
(som)
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

