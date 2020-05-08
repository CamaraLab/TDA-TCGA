#' Summarizes and identifies significant mutations from filter_genes
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, and nerve complexes
#' @param q_threshold_localization threshold q value for gene localization using Benjamini-Hochberg procedure to control false discovery rate
#'
#' @return Returns significant nonsynonymously mutated genes
#'
#' @export

identify_significant_mutations <- function(TDAmut_object, q_threshold_localization = 0.15) {

  if (is_empty(TDAmut_object@gene_scores)){
    stop('Please run compute_complexes, compute_mut_load, and filter_genes first to populate object with nerve complexes, formatted mutation data, and gene localization scores')
  }

  gene_scores_list <- TDAmut_object@gene_scores
  neg_cor_list <- TDAmut_object@negative_correlations
  num_complexes <- length(gene_scores_list)
  filtered_genes <- TDAmut_object@filtered_genes

  ######## FILTERING BY MEDIAN Q VALUE THRESHOLD ########
  
  # Rearranging scores by gene instead of by nerve complex
  collapsed_cor <- bind_rows(neg_cor_list)
  genes_cor <- unique(collapsed_cor$Gene)
  neg_cor_q_bygene <- matrix(collapsed_cor$q, length(genes_cor), num_complexes, dimnames = list(genes_cor, NULL), byrow = FALSE) %>% as.data.frame # genes x complexes
  
  collapsed_loc <- bind_rows(gene_scores_list)
  genes_loc <- unique(collapsed_loc$Gene)
  localization_q_bygene <- matrix(collapsed_loc$q0, length(genes_loc), num_complexes, dimnames = list(genes_loc, NULL), byrow = FALSE) %>% as.data.frame
  
  neg_cor_genes <- filtered_genes[!(filtered_genes %in% rownames(localization_q_bygene))] # Genes displaying negative correlations
  
  sig_genes <- localization_q_bygene[apply(localization_q_bygene, 1, function(x) any(median(x) < q_threshold_localization)), ] # Significantly localized genes
  
  ######## SUMMARIZING RESULTS ########
  
  # Genes with significant localization
  message('\nThe following mutated genes were identified as significant given the input localization q value threshold: ')
  
  message(paste("'", rownames(sig_genes), "'", collapse = ", ", sep = ""))
  
  # Genes with negative correlation
  message('\nThe following mutated genes were discarded due to negative correlations in expression and mutation profiles: ')
  
  message(paste("'", neg_cor_genes, "'", collapse = ", ", sep = ""))
  
  # Genes with nonsignificant localization
  excluded_genes <- localization_q_bygene[!(rownames(localization_q_bygene) %in% rownames(sig_genes) )]
  
  message('\nThe following mutated genes were discarded due to localization since their median q value below the input threshold: ')
  
  message(paste("'", rownames(excluded_genes), "'", collapse = ", ", sep = ""))

  
  TDAmut_object@significant_genes <- sig_genes
  
  return(TDAmut_object)
}
