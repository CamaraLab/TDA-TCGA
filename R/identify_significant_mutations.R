#' Identifies significant genes from asses_mutations
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, and nerve complexes.
#' @param top_nonsyn_fraction number of genes to keep with greatest nonsyn/nonsyn+syn fraction. By default is 350.
#' @param q_threshold threshold q value using Benjamini-Hochberg procedure to control false discovery rate
#'
#' @return Returns table of significant nonsynonymously mutated genes
#'
#' @export

identify_significant_mutations <- function(TDAmut_object, top_nonsyn_fraction = 350, q_threshold = 0.15) {

  if (is_empty(TDAmut_object@nerve_complexes) || is_empty(TDAmut_object@nonsyn_mutations) || is_empty(TDAmut_object@gene_scores)){
    stop('Please run compute_complexes, compute_mut_load, and assess_mutations first to populate object with nerve complexes, formatted mutation data, and gene localization scores')
  }
  
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  syn_mat <- TDAmut_object@syn_mutations
  gene_scores_table <- TDAmut_object@gene_scores

  # Apply thresholds to p value, q value, and nonsyn fraction
  nonsyn_frac <- colSums(nonsyn_mat) / colSums(nonsyn_mat + syn_mat[ , colnames(nonsyn_mat)])
  top_frac <- colnames(nonsyn_mat[ , order(nonsyn_frac, decreasing = TRUE)][ , 1:top_nonsyn_fraction])
  
  # gene_scores_table <- gene_scores_table[gene_scores_table$Gene %in% top_frac, ]
  # gene_scores_table <- gene_scores_table[gene_scores_table$p-value <= 0.05,] 
  # gene_scores_table <- gene_scores_table[gene_scores_table$q-value < q_threshold, ]

  sig_genes <- lapply(gene_scores_table, function(x){
    x <- x[x$p0 <= 0.05, ]
    x <- x[x$q0 < q_threshold, ]
    x <- x[x$Gene %in% top_frac, ]
    return(x)
  })
  
  TDAmut_object@significant_genes <- sig_genes
  
  return(TDAmut_object)
}
