library(RayleighSelection)
library(dplyr)

compute_mut_load <- function(exp_table, mut_table, mut_objects, hypermut_cutoff=FALSE, nonsyn_muts=NA, syn_muts=NA, graph=FALSE, num_cores=1) {

  exp_table <- read.csv("/home/rstudio/documents/Messy_test_data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F)
  mut_table <- read.csv('/home/rstudio/documents/TDA-TCGA/Test_Data/LGG_Muts.txt',row.names=1,header=T,stringsAsFactors = F)
  # exp_table <- read.csv(exp_table, row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?"))
  # mut_table <- read.csv(mut_table, row.names=1, header=T, stringsAsFactors = F, na.strings=c("NA","NaN", " ", "?")
                    
  samples <- sort(unique(mut_table$Sample))
  gene_names <- sort(unique(mut_table$Gene))
  
  # Assigning syn and nonsyn mutations
  if(is.na(nonsyn_muts) || is.na(syn_muts)){
    nonsyn_type <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site",      
                     "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins",
                     "In_Frame_Ins", "Nonstop_Mutation")
    
    nonsyn_muts <- mut_table[mut_table$Type %in% nonsyn_type,]
    syn_muts <- mut_table[!(mut_table$Type %in% nonsyn_type),]
  }
  
  # Grouping by gene and sample names
  nonsyn_mat <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
  dimnames(nonsyn_mat) <- list(samples,gene_names)
  syn_mat <- nonsyn_mat
  
  t_nonsyn <- with(nonsyn_muts,table(Sample,Gene))
  nonsyn_mat[rownames(t_nonsyn),colnames(t_nonsyn)] <- t_nonsyn
  t_syn <- with(syn_muts,table(Sample,Gene))
  syn_mat[rownames(t_syn),colnames(t_syn)] <- t_syn
  
  # Harmonizing names between expression and mutation data
  if(!all(rownames(exp_table) %in% rownames(nonsyn_mat)) ) {
    missing_samples <- rownames(exp_table[!(rownames(exp_table) %in% rownames(nonsyn_mat)),])
    message(paste0("Excluding following samples from nonsynonymous mutation table: ", missing_samples))
    
    # newmat <- matrix(0,length(missing_samples),ncol(nonsyn_mat))
    # dimnames(newmat) <- list(missing_samples,colnames(nonsyn_mat))
    # nonsyn_mat <- rbind(nonsyn_mat,newmat)
    # nonsyn_mat <- nonsyn_mat[order(rownames(nonsyn_mat)),]
  }
  
  if(!all(rownames(exp_table) %in% rownames(syn_mat)) ) {
    missing_samples <- rownames(exp_table[!(rownames(exp_table) %in% rownames(syn_mat)),])
    message(paste0("Excluding following samples from synonymous mutation table: ", missing_samples))
    
    # newmat <- matrix(0,length(missing_samples),ncol(syn_mat))
    # dimnames(newmat) <- list(missing_samples,colnames(syn_mat))
    # syn_mat <- rbind(syn_mat,newmat)
    # syn_mat <- syn_mat[order(rownames(syn_mat)),]
  }
  
  # Creating binary nonsynonymous mutation matrix
  #  nonsyn_bin<-ifelse(nonsyn_mat>0,1,0)

  # Computing total mutational load
  mutload <- rowSums(syn_mat+nonsyn_mat)
  
  # Subsampling mutations
  if(hypermut_cutoff != FALSE){
    
    rescale_boundary <- 10^hypermut_cutoff
    above_cutoff <- mutload[mutload>rescale_boundary]
    below_cutoff <- mutload[mutload<=rescale_boundary]
    subsample <- floor(median(above_cutoff) / median(below_cutoff))
    
    hypermut_samples <- names(above_cutoff)
    muts <- mut_table[mut_table$Sample %in% hypermut_samples,]
    to_keep <- floor(nrow(muts) / subsample)
    swapped <- sample(nrow(muts),to_keep)
    muts_sub <- muts[swapped,]
    mut_table <- mut_table[!(mut_table$Sample %in% hypermut_samples),]
    mut_table <- rbind(mut_table,muts_sub)
    mut_table <- mut_table[order(mut_table$Sample),]
    
    samples <- sort(unique(mut_table$Sample))
    gene_names <- sort(unique(mut_table$Gene))
    
    nonsyn_muts <- mut_table[mut_table$Type %in% nonsyn_type,]
    syn_muts <- mut_table[!(mut_table$Type %in% nonsyn_type),]
    
    nonsyn_mat <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
    dimnames(nonsyn_mat) <- list(samples,gene_names)
    syn_mat <- nonsyn_mat
    
    t_nonsyn <- with(nonsyn_muts,table(Sample,Gene))
    nonsyn_mat[rownames(t_nonsyn),colnames(t_nonsyn)] <- t_nonsyn
    nonsyn_mat <- nonsyn_mat[order(rownames(nonsyn_mat)), order(colnames(nonsyn_mat))]
    t_syn <- with(syn_muts,table(Sample,Gene))
    syn_mat[rownames(t_syn),colnames(t_syn)] <- t_syn
    syn_mat <- syn_mat[order(rownames(syn_mat)), order(colnames(syn_mat))]
    
    mutload <- rowSums(syn_mat+nonsyn_mat)
  }
  
  # Computing localization of mutational load
  mutload_mat <- t(mutload) %>% as.data.frame
  num <- length(mut_objects)
  pvals <- vector('list',num)
  
  for(i in seq(1,num,1)){
    pvals[[i]] <- rayleigh_selection(mut_objects[[i]], mutload_mat, num_perms = 2000, num_cores=13,one_forms=FALSE)
  }
  
  if(graph=TRUE){
    
    # Total mutational load histogram
    hist(log10(mutload),breaks=100,plot=TRUE,main="LGG Mutational Load")
    
    # Localization across grid of Mapper parameters
    intervals = seq(10, 80, by = 10) 
    #percents = seq(15,85,by = 10)
    percents <- c(33,60,71.4,77.8,81.8,84.6,86.7,88.2)
    pval_data <- matrix(sapply(pvals_sub,'[',2),8,8,byrow=TRUE) %>% as.data.frame
    colnames(pval_data) <- percents
    dat2 <- cbind(intervals,pval_data) %>% as.matrix
    mode(dat2) = 'numeric'
    dat2 <- dat2 %>% as.data.frame
    data <- melt(dat2,id.vars='intervals',measure.vars = c(as.character(percents)))
    data <- data.frame(Intervals = rep(intervals,each=8), Percents = rep(percents,8), Values=as.numeric(data$value))
    
    ggplot(data, aes(x=factor(Intervals),y=factor(Percents),fill=Values,label=round(data$Values,3))) +
      geom_tile(alpha=0.7) + theme_minimal()+ geom_text(size=3) +
      scale_fill_gradient2(low='red',high='blue',mid='purple',midpoint=0.5,guide_legend(title="p-value")) +
      ggtitle("Localization across Mapper Complexes") +
      xlab("Number of Intervals") + ylab ("Percent Overlap") +
      theme(axis.text=element_text(size= 8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold"),panel.grid.major = element_blank()) +
      theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) +
      scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      coord_equal()
  }
  
  mutload_object <- vector('list',length = 3)
  mutload_object[[1]] <- nonsyn_mat
  mutload_object[[2]] <- syn_mat
  mutload_object[[3]] <- mutload
  class(mutload_object) <- 'Mutational_Load'
  return(mutload_object)
  
  
  
  ###############################################################################
  
  # Subsampling
  if(!is.null(hypermut_cutoff)){
    
    rescale_boundary <- 10^hypermut_cutoff
    above_cutoff <- mutload[mutload>rescale_boundary]
    below_cutoff <- mutload[mutload<=rescale_boundary]
    subsample <- floor(median(above_cutoff) / median(below_cutoff))
    
    hypermut_samples <- names(above_cutoff)
    muts <- mut_table[mut_table$Sample %in% hypermut_samples,]
    to_keep <- floor(nrow(muts) / subsample)
    swapped <- sample(nrow(muts),to_keep)
    muts_sub <- muts[swapped,]
    mut_table <- mut_table[!(mut_table$Sample %in% hypermut_samples),]
    mut_table <- rbind(mut_table,muts_sub)
    mut_table <- mut_table[order(mut_table$Sample),]

    samples <- sort(unique(mut_table$Sample))
    gene_names <- sort(unique(mut_table$Gene))
    
    nonsyn_muts <- mut_table[mut_table$Type %in% nonsyn_type,]
    syn_muts <- mut_table[!(mut_table$Type %in% nonsyn_type),]
    
    nonsyn_mat <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
    dimnames(nonsyn_mat) <- list(samples,gene_names)
    syn_mat <- nonsyn_mat
    
    t_nonsyn <- with(nonsyn_muts,table(Sample,Gene))
    nonsyn_mat[rownames(t_nonsyn),colnames(t_nonsyn)] <- t_nonsyn
    nonsyn_mat <- nonsyn_mat[order(rownames(nonsyn_mat)), order(colnames(nonsyn_mat))]
    t_syn <- with(syn_muts,table(Sample,Gene))
    syn_mat[rownames(t_syn),colnames(t_syn)] <- t_syn
    syn_mat <- syn_mat[order(rownames(syn_mat)), order(colnames(syn_mat))]
    
    if(!all(rownames(exp_table) %in% rownames(nonsyn_mat)) ) {
      missing_samples <- rownames(exp_table[!(rownames(exp_table) %in% rownames(nonsyn_mat)),])
      newmat <- matrix(0,length(missing_samples),ncol(nonsyn_mat))
      dimnames(newmat) <- list(missing_samples,colnames(nonsyn_mat))
      nonsyn_mat <- rbind(nonsyn_mat,newmat)
      nonsyn_mat <- nonsyn_mat[order(rownames(nonsyn_mat)),]
    }
    
    if(!all(rownames(exp_table) %in% rownames(syn_mat)) ) {
      missing_samples <- rownames(exp_table[!(rownames(exp_table) %in% rownames(syn_mat)),])
      newmat <- matrix(0,length(missing_samples),ncol(syn_mat))
      dimnames(newmat) <- list(missing_samples,colnames(syn_mat))
      syn_mat <- rbind(syn_mat,newmat)
      syn_mat <- syn_mat[order(rownames(syn_mat)),]
    }
    
    # Plotting histogram of mutational load
    sub_mutload <- rowSums(syn_mat+nonsyn_mat)
    hist(log10(sub_mutload),breaks=100,plot=TRUE,main="Subsampled LGG Mutational Load")
    sub_mutload_mat <- t(sub_mutload) %>% as.data.frame
    
    # Plotting localization of mutational load
    num <- length(mut_objects)
    pvals_sub <- vector('list',num)
    
    for(i in seq(1,num,1)){
      pvals_sub[[i]] <- rayleigh_selection(mut_objects[[i]],sub_mutload_mat,num_cores=13,one_forms=FALSE, num_perms = 2000 )
    }
  }

}
