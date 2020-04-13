#' Compute each gene's Laplacian score, p value, and q value via BH procedure for controlling false discovery rate across grid of nerve complexes
#' 
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes
#' @param top_nonsyn_fraction number of genes to keep with greatest nonsyn/nonsyn+syn fraction. By default is 350.
#' @param num_permutations number of permutations used to create null distribution. By default is 10000
#' @param num_cores number of cores to be used in computation. By default is 1.
#' 
#' @return Returns TDAmut object populated with genes' Laplacian scores, p values, and q values
#' 
#' @export

assess_mutations <- function(TDAmut_object, top_nonsyn_fraction = 350, num_permutations = 10000, num_cores = 1){
  
  nerve_complexes <- TDAmut_object@nerve_complexes
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  
  # Nonsynonymous fraction threshold
  nonsyn_frac <- colSums(nonsyn_mat) / colSums(nonsyn_mat + syn_mat[ , colnames(nonsyn_mat)])
  nonsyn_mat <- nonsyn_mat[ , nonsyn_frac > top_nonsyn_fraction]
  
  nonsyn_features <- t(nonsyn_mat) %>% as.data.frame
  num <- length(nerve_complexes)
  pvals <- vector('list',num)
  
  for(i in seq(1,num,1)){
    pvals[[i]] <- rayleigh_selection(nerve_complexes[[i]], nonsyn_features, num_perm = num_permutations, num_cores = num_cores, one_forms = FALSE)
  }
  
  TDAmut_object@gene_localization <- pvals
  
}