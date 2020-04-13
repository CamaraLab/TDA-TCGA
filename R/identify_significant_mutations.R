#' Identifies significant genes from asses_mutations
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, and nerve complexes.
#' @param freq_threshold threshold mutation frequency. Genes below this value are not considered in analysis. By default is 0.02
#' @param q_threshold threshold q value using Benjamini-Hochberg procedure to control false discovery rate
#' @param num_cores number of cores to be used in computation. By default is 1.
#'
#' @return Returns table of significant nonsynonymously mutated genes
#'
#' @export

identify_significant_mutations <- function(TDAmut_object, freq_threshold = 0.02, q_threshold = 0.15, num_cores = 1) {

  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  pvals <- TDAmut_object@gene_localization

  gene_local_table <- data.frame("Gene"= rownames(nonsyn_features), "Laplacian Score" = sapply(pvals,'[',1), "p-value" = sapply(pvals,'[',2), "q-value" = sapply(pvals,'[',3))
  
  # Apply thresholds to p value, q value, and frequency
  gene_local_table <- gene_local_table[gene_local_table$p-value <= 0.05,]
  gene_local_table <- gene_local_table[gene_local_table$q-value < q_threshold, ]
  
  nonsyn_bin<-ifelse(nonsyn_mat>0,1,0) %>% as.data.frame
  freq <- colMeans(nonsyn_bin) %>% as.data.frame
  freq_genes <- colnames(nonsyn_mat[ , freq > freq_threshold])
  gene_local_table <- gene_local_table[gene_local_table$Gene %in% freq_genes, ]
  
  TDAmut_object@significant_genes <- gene_local_table
  
  return(gene_local_table)
}
