#' Assesses localization of nonsynonymous mutations across grid of nerve complexes
#'
#' @param mut_objects object output by create_mut_object
#' @param mutload_object object output by compute_mut_load
#' @param freq_threshold threshold mutation frequency. Genes below this value are not considered in analysis. By default is 0.02
#' @param top_nonsyn_fraction number of genes to keep with greatest nonsyn/nonsyn+syn fraction. By default is 350.
#' @param permutations number of permutations used to create null distribution. By default is 10000
#' @param q_threshold threshold q value using Benjamini-Hochberg procedure to control false discovery rate
#' @param num_cores number of cores to be used in computation. By default is 1.
#'
#' @return Returns table of significant nonsynonymously mutated genes
#'
#' @export

identify_drivers <- function(mut_objects, mutload_object, freq_threshold = 0.02, top_nonsyn_fraction = 350, permutations = 5000, q_threshold = 0.15, num_cores=1) {

  nonsyn_mat <- mutload_objects[[1]]

  samples <- rownames(nonsyn_mat)
  genes <- colnames(nonsyn_mat)

  # Freq threshold
  nonsyn_bin<-ifelse(nonsyn_mat>0,1,0) %>% as.data.frame
  freq <- colMeans(nonsyn_bin) %>% as.data.frame
  nonsyn_mat <- nonsyn_mat[,freq > freq_threshold]

  # NonSynonymous fraction threshold
  nonsyn_frac <- colSums(nonsyn_mat) / colSums(nonsyn_mat + syn_mat[,colnames(nonsyn_mat)])
  nonsyn_mat <- nonsyn_mat[,nonsyn_frac > top_nonsyn_fraction]

  # Assess localization of genes passing threshold
  nonsyn_features <- t(nonsyn_mat) %>% as.data.frame
  num <- length(mut_objects)
  pvals <- vector('list',num)

  for(i in seq(1,num,1)){
    pvals[[i]] <- rayleigh_selection(mut_objects[[i]], nonsyn_features, num_perms = permutations, num_cores=13,one_forms=FALSE)
  }

  # Keep genes with p < 0.05 and q < 0.15
  gene_local_table <- data.frame("Gene"= rownames(nonsyn_features), "Laplacian Score" = sapply(pvals,'[',1), "p-value" = sapply(pvals,'[',2), "q-value" = sapply(pvals,'[',3))

  gene_local_table <- gene_local_table[gene_local_table$p-value < 0.05,]
  gene_local_table <- gene_local_table[gene_local_table$q-value < q_threshold, ]
}
