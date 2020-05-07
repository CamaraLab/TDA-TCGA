#' Computes mutational load from input mutation table and assesses localization across grid of nerve complexes
#'
#' @param TDAmut_object Object of class TDAmut with expression data, mutation data, and nerve complexes
#' @param hypermut_cutoff threshold in log10 scale for downsampling mutations. By default is FALSE.
#' @param min_mutload minimum threshold in log10 scale for mutational load of samples. Samples below threshold are saved and to be used to recompute nerve complexes in compute_complexes
#' @param nonsyn_muts optional input of nonsynonymous mutations in samples. By default is NA.
#' @param syn_muts optional input of synonymous mutations in samples. By default is NA.
#' @param num_permutations number of permutations used to create null distribution. By default is 2000
#' @param num_cores number of cores to be used in computation. By default is 1.
#' @param seed integer specifying the seed used to initialize the generator of permutations. By default is 121

#'
#' @return Appends object of class TDAmut with list of pvalues indicating localization of mutational load and 
#' data frames of nonsynonymous mutations, synonymous mutations, and mutational load of samples
#'
#' @export

# library(RayleighSelection)
# library(dplyr)

compute_mut_load <- function(TDAmut_object, hypermut_cutoff = FALSE, min_mutload = NULL, nonsyn_muts = NULL, syn_muts = NULL, num_permutations = 1000, num_cores = 1, seed = 121) {
  
  if (is_empty(TDAmut_object@nerve_complexes) || is_empty(TDAmut_object@mutation_table)){
    stop('Run create_TDAmut_object and compute_complexes first to populate object with mutation data and nerve_complexes')
  }
  
  ######## CONSOLIDATE MUTATION TABLE INTO DATA FRAMES ########
  
  split_mut_data <- function(TDAmut_object) {
    
    mut_table <- TDAmut_object@mutation_table
    
    samples <- sort(unique(mut_table$Sample))
    gene_names <- sort(unique(mut_table$Gene))
    
    # Reformatting table into matrix
    mut_mat <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
    dimnames(mut_mat) <- list(samples,gene_names)
    t_mut <- with(mut_table,table(Sample,Gene))
    mut_mat[rownames(t_mut),colnames(t_mut)] <- t_mut
    mut_mat <- mut_mat[order(rownames(mut_mat)), order(colnames(mut_mat))]
    
    # Assigning syn and nonsyn mutations
    if(is.null(nonsyn_muts) || is.null(syn_muts)){
      nonsyn_type <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
                       "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins",
                       "In_Frame_Ins", "Nonstop_Mutation")
      
      nonsyn_muts <- mut_table[mut_table$Type %in% nonsyn_type,]
      syn_muts <- mut_table[!(mut_table$Type %in% nonsyn_type),]
    }
    
    # Reformatting nonsynonymous and synonymous mutations into matrices
    nonsyn_mat <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
    dimnames(nonsyn_mat) <- list(samples,gene_names)
    syn_mat <- nonsyn_mat
    
    t_nonsyn <- with(nonsyn_muts,table(Sample,Gene))
    nonsyn_mat[rownames(t_nonsyn),colnames(t_nonsyn)] <- t_nonsyn
    nonsyn_mat <- nonsyn_mat[order(rownames(nonsyn_mat)), order(colnames(nonsyn_mat))]
    t_syn <- with(syn_muts,table(Sample,Gene))
    syn_mat[rownames(t_syn),colnames(t_syn)] <- t_syn
    syn_mat <- syn_mat[order(rownames(syn_mat)), order(colnames(syn_mat))]
    
    TDAmut_object@mutation_matrix <- mut_mat
    TDAmut_object@nonsyn_mutations <- nonsyn_mat
    TDAmut_object@syn_mutations <- syn_mat
    TDAmut_object@mutational_load <- rowSums(syn_mat + nonsyn_mat)
    
    return(TDAmut_object)
  }
  
  TDAmut_object <- split_mut_data(TDAmut_object)

  ######## SUBSAMPLING MUTATIONS ########
  
  if(hypermut_cutoff != FALSE){
    mutload <- TDAmut_object@mutational_load
    mut_table <- TDAmut_object@mutation_table
    
    # Equalizing median mutational load between hypermutated samples and other samples
    rescale_boundary <- 10^hypermut_cutoff
    above_cutoff <- mutload[mutload>rescale_boundary]
    below_cutoff <- mutload[mutload<=rescale_boundary]
    subsample <- floor(median(above_cutoff) / median(below_cutoff))

    hypermut_samples <- names(above_cutoff)
    message('The following samples surpass the mutational load threshold of ', 10^hypermut_cutoff, ': ', paste("'", hypermut_samples, "'", collapse = ", ", sep = ""))
    muts <- mut_table[mut_table$Sample %in% hypermut_samples,]
    to_keep <- floor(nrow(muts) / subsample)
    swapped <- sample(nrow(muts),to_keep)
    muts_sub <- muts[swapped,]
    mut_table <- mut_table[!(mut_table$Sample %in% hypermut_samples),]
    mut_table <- rbind(mut_table,muts_sub)
    mut_table <- mut_table[order(mut_table$Sample),]
    
    TDAmut_object@mutation_table <- mut_table
    
    TDAmut_object <- split_mut_data(TDAmut_object)
  }
  
  ######## COMPUTING LOCALIZATION OF MUTATIONAL LOAD ########
  
  nerve_complexes <- TDAmut_object@nerve_complexes
  mutload <- TDAmut_object@mutational_load
  
  if (min_mutload != FALSE){
    min_samples <- names(mutload[mutload < 10^min_mutload])
    TDAmut_object@min_mutated_samples <- min_samples
    message('The following samples have a mutational load less than ', floor(10^min_mutload), ": ", paste("'", min_samples, "'", collapse = ", ", sep = ""))
    message('Consider using compute_complexes to recompute nerve complexes without these samples')
  }
  
  mutload_mat <- t(mutload) %>% as.data.frame
  num_complexes <- length(nerve_complexes)
  pvals <- vector('list', num_complexes)

  for(i in seq(1, num_complexes, 1)){
    message(paste("Assessing localization in nerve complex", i, "of", num_complexes))
    pvals[[i]] <- rayleigh_selection(nerve_complexes[[i]], mutload_mat, num_perm = num_permutations, num_cores = num_cores, one_forms = FALSE, seed = seed)
  }
  
  TDAmut_object@mutational_load_localization <- pvals
  
  return(TDAmut_object)
  
}