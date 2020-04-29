#' Summarizes and identifies significant mutations from filter_genes
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, and nerve complexes
#' @param q_threshold_localization threshold q value for gene localization using Benjamini-Hochberg procedure to control false discovery rate
#' @param q_threshold_negative_correlation threshold q value for negative correlations using Benjamini-Hochberg procedure to control false discovery rate. Genes with q above this threshold display negative correlations between expression and mutation data.
#' @param remove_negative_correlation removes genes with significant negative correlations between expression and mutation data. Must have negative_correlations = TRUE in filter_genes. By default is TRUE.
#'
#' @return Returns table of significant nonsynonymously mutated genes
#'
#' @export

identify_significant_mutations <- function(TDAmut_object, q_threshold_localization = 0.15, remove_negative_correlation = TRUE, q_threshold_negative_correlation = 0.80) {

  if (is_empty(TDAmut_object@gene_scores)){
    stop('Please run compute_complexes, compute_mut_load, and filter_genes first to populate object with nerve complexes, formatted mutation data, and gene localization scores')
  }
  
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  syn_mat <- TDAmut_object@syn_mutations
  gene_scores_list <- TDAmut_object@gene_scores

  ######## FILTERING BY p VALUE, q VALUE, and NEGATIVE CORRELATIONS ########

  # gene_scores_table <- gene_scores_table[gene_scores_table$Gene %in% top_frac, ]
  # gene_scores_table <- gene_scores_table[gene_scores_table$p-value <= 0.05,] 
  # gene_scores_table <- gene_scores_table[gene_scores_table$q-value < q_threshold, ]

  sig_genes <- lapply(gene_scores_list, function(x){
    x <- x[x$p0 <= 0.05, ]
    x <- x[x$q0 < q_threshold, ]
    return(x)
  })
  
  if (remove_negative_correlation == TRUE){
    neg_cor_list <- TDAmut_object@negative_correlations 
    sig_genes_no_cor <- sig_genes
    
    neg_cor_genes <- lapply(neg_cor_list, function(x){
      x <- x$Gene[x$q > q_threshold_negative_correlation]
    }
    
    for (i in length(nerve_complexes)){
      genetable <- sig_genes[[i]]
      cortable <- neg_cor_list[[i]]
      sig_genes_no_cor[[i]] <- genetable[!(genetable$Gene %in% cortable$Gene), ]
    }
    
  }

  ######## SUMMARIZING SIGNIFICANT GENES ########
  
  #summary of genes removed by q value and p value
  
  #summary of genes removed by neg_corr (look at sig_genes_no_cor)
  
  
  TDAmut_object@significant_genes <- sig_genes
  
  return(TDAmut_object)
}
