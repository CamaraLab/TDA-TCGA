library(RayleighSelection)
library(dplyr)

compute_mut_load <- function(exp_table, mut_table, mut_objects, hypermut_cutoff, nonsyn_muts=NA, syn_muts=NA, num_cores=1) {

  exp_table <- (read.csv(exp_table, row.names=1, header=T, stringsAsFactors=F))
  mut_table <- (read.csv(mut_table, row.names=1, header=T, stringsAsFactors=F))
                         
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
  t_syn <- with(syn_muts,table(Sample,Gene,Type))
  syn_mat[rownames(t_syn),colnames(t_syn)] <- t_syn
  
  # Harmonizing names between expression and mutation data
  if(!all(rownames(exp_table) %in% rownames(nonsyn_mat)) ) {
    missing_samples <- rownames(exp_table[!(rownames(exp_table) %in% rownames(nonsyn_mat)),])
    newmat <- matrix(0,length(missing_samples),ncol(nonsyn_mat))
    colnames(newmat) <- colnames(nonsyn_mat)
    nonsyn_mat <- rbind(nonsyn_mat,newmat)
    nonsyn_mat <- nonsyn_mat[order(rownames(nonsyn_mat)),]
  }
  
  if(!all(rownames(exp_table) %in% rownames(syn_mat)) ) {
    missing_samples <- rownames(exp_table[!(rownames(exp_table) %in% rownames(syn_mat)),])
    newmat <- matrix(0,length(missing_samples),ncol(syn_mat))
    colnames(newmat) <- colnames(syn_mat)
    syn_mat <- rbind(syn_mat,newmat)
    syn_mat <- nonsyn_mat[order(rownames(syn_mat)),]
  }
  
  # Creating binary nonsynonymous mutation matrix for later use
  nonsyn_bin<-ifelse(nonsyn_mat>0,1,0)

  # Computing and plotting total mutational load
  mutload <- rowSums(syn_mat+nonsyn_mat)
  hist(log10(mutload),breaks=100,plot=TRUE,main="LGG Mutational Load")
  
  # Plotting localization of mutational load
  mutload_mat <- t(mutload) %>% as.data.frame
  num <- length(mut_objects)
  pvals <- vector('list',num)
  
  for(i in seq(1,num,1)){
    pvals[[i]] <- rayleigh_selection(mut_objects[[i]],mutload_mat,num_cores=10,one_forms=FALSE)
  }
  
  # Subsampling mutations
  if(!is.null(hypermut_cutoff)){
    
    rescale_boundary <- 10^hypermut_cutoff
    above_cutoff <- mutload[mutload>rescale_boundary]
    below_cutoff <- mutload[mutload<=rescale_boundary]
    subsample <- floor(median(above_cutoff) / median(below_cutoff))
    
    hypermut_samples <- names(above_cutoff)
    muts <- mut_table[mut_table$Sample %in% hypermut_samples,]
    swapped <- sample(nrow(muts))
    muts <- muts[swapped,]
    index_hold <- seq.int(1,nrow(muts),by=subsample)
    muts_sub <- muts[index_hold,]
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
    nonsyn_mat <- nonsyn_mat[sort(rownames(nonsyn_mat)), sort(colnames(nonsyn_mat))]
    t_syn <- with(syn_muts,table(Sample,Gene,Type))
    syn_mat[rownames(t_syn),colnames(t_syn)] <- t_syn
    syn_mat <- syn_mat[sort(rownames(syn_mat)), sort(colnames(syn_mat))]
    
    sub_mutload <- rowSums(syn_mat+nonsyn_mat)
    hist(log10(sub_mutload),breaks=100,plot=TRUE,main="Subsampled LGG Mutational Load")
    
    if(! all (rownames(exp_table) %in% names(sub_mutload)) ) {
      missing_samples <- rownames(exp_table[!(rownames(exp_table) %in% names(sub_mutload)),])
      #missing_samples <- c("TCGA-DU-7014-01A","TCGA-DU-A7TI-01A","TCGA-HW-7493-01A","TCGA-TQ-A7RK-02B")
      newmat <- matrix(0,1,length(missing_samples))
      colnames(newmat) <- missing_samples
      sub_mutload_mat <- cbind(t(sub_mutload) %>% as.data.frame, newmat)
      sub_mutload_mat <- sub_mutload_mat[,order(colnames(sub_mutload_mat))] %>% as.data.frame
    }
      else{
        sub_mutload_mat <- t(mutload) %>% as.data.frame
      }
    
    # Plotting localization of mutational load
    num <- length(mut_objects)
    pvals_sub <- vector('list',num)
    
    for(i in seq(1,num,1)){
      pvals_sub[[i]] <- rayleigh_selection(mut_objects[[i]],sub_mutload_mat,num_cores=12,one_forms=FALSE)
    }
    
    
    
  }

}
