#' Compute each gene's Laplacian score, p value, and q value via BH procedure for controlling false discovery rate across grid of nerve complexes
#' 
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes
#' @param freq_threshold threshold mutation frequency. Genes below this value are not considered in analysis. By default is 0.02
#' @param num_permutations number of permutations used to create null distribution. By default is 10000
#' @param num_cores number of cores to be used in computation. By default is 1.
#' 
#' @return Returns TDAmut object populated with genes' Laplacian scores, p values, and q values
#' 
#' @export

assess_mutations <- function(TDAmut_object, freq_threshold = 0.01, num_permutations = 10000, num_cores = 5){
  
  if (is_empty(TDAmut_object@nerve_complexes) || is_empty(TDAmut_object@nonsyn_mutations)){
    stop('Run compute_complexes and compute_mut_load first to populate object with nerve complexes and formatted mutation data')
  }
  
  nerve_complexes <- TDAmut_object@nerve_complexes
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  
  # Freq Threshold
  nonsyn_bin <- ifelse(nonsyn_mat > 0, 1, 0) %>% as.data.frame
  freq <- colMeans(nonsyn_bin) %>% as.data.frame
  nonsyn_mat <- nonsyn_mat[ , freq > freq_threshold]
  
  nonsyn_features <- t(nonsyn_mat) %>% as.data.frame
  num_complexes <- length(nerve_complexes)
  gene_scores <- vector('list', num_complexes)
  
  for(i in seq(1, num_complexes, 1)){
    message(paste("Assessing genes in nerve complex", i, "of", num_complexes))
    gene_scores[[i]] <- rayleigh_selection(nerve_complexes[[i]], nonsyn_features, num_perm = num_permutations, num_cores = num_cores, one_forms = FALSE)
    gene_scores[[i]] <- data.frame("Gene"= rownames(nonsyn_features), "R" = sapply(gene_scores[i], '[', 1), "p" = sapply(gene_scores[i], '[', 2), "q" = sapply(gene_scores[i], '[', 3))
  }
  
  TDAmut_object@gene_scores <- gene_scores
  
  return(TDAmut_object)
  
}