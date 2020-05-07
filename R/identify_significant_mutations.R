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
  num_complexes <- length(gene_scores_list)

  ######## FILTERING BY p, q, NEG CORRELATIONS (PER COMPLEX) ########
# 
#   # gene_scores_table <- gene_scores_table[gene_scores_table$Gene %in% top_frac, ]
#   # gene_scores_table <- gene_scores_table[gene_scores_table$p-value <= 0.05,] 
#   # gene_scores_table <- gene_scores_table[gene_scores_table$q-value < q_threshold, ]
#   
#   sig_genes <- lapply(gene_scores_list, function(x){
#     # x <- x[x$p0 <= 0.05, ]
#     excluded_localization <- x$Gene[x$q0 > q_threshold] # maybe change
#     x <- x[x$q0 <= q_threshold, ]
#     return(x)
#   })
#   
#   if (remove_negative_correlation == TRUE){
#     neg_cor_list <- TDAmut_object@negative_correlations 
#     sig_genes_no_cor <- sig_genes
#     
#     excluded_neg_cor <- lapply(neg_cor_list, function(x){
#       # x <- x$Genes[x$q > q_threshold_negative_correlation]
#       x <- x$genes.gene[x$q > q_threshold_negative_correlation]
#     })
#     
#     for (i in 1:num_complexes{
#       genetable <- sig_genes[[i]]
#       excluded_neg_cor_genes <- excluded_neg_cor[[i]]
#       sig_genes_no_cor[[i]] <- genetable[!(genetable$Gene %in% excluded_neg_cor_genes), ]
#     }
#     
#   }
  
  ######## FILTERING NEGATIVE CORRELATION (FROM PAPER) ########
  
  # Rearranging scores by gene instead of by nerve complex
  collapsed <- bind_rows(neg_cor_list)
  neg_cor_q_bygene <- matrix(collapsed$q, 6, 36, dimnames = list(unique(collapsed$genes.gene., NULL)), byrow = FALSE) %>% as.data.frame # genes x complexes
  #q_bygene <- matrix(collapsed$q, 36, 6, dimnames = list(NULL, unique(collapsed$genes.gene.)), byrow = TRUE) %>% as.data.frame %>% as.list # genes as elements of list
  #q_bygene <- matrix(collapsed$q, 36, 6, dimnames = list(NULL, unique(collapsed$genes.gene.)), byrow = TRUE) %>% as.data.frame # complexes x genes
  
  collapsed <- bind_rows(gene_scores_list)
  localization_q_bygene <- matrix(collapsed$q0, 6, 36, dimnames = list(unique(collapsed$gene., NULL)), byrow = FALSE) %>% as.data.frame
  
  neg_cor_genes <- neg_cor_q_bygene[apply(neg_cor_q_bygene, 1, function(x) any(median(x) > q_threshold_negative_correlation)), ]
  
  genes_no_cor <- localization_q_bygene[!(rownames(localization_q_bygene) %in% rownames(neg_cor_genes)), ] # Not considering genes with neg correlation q val above threshold
  
  sig_genes <- genes_no_cor[apply(genes_no_cor, 1, function(x) any(median(x) < q_threshold_localization)), ]
  
  ######## SUMMARIZING SIGNIFICANT GENES ########
  
  # Summary of significant genes
  message('\nThe following mutated genes were identified as significant given the input q value threshold: ')
  
  message(paste("'", rownames(sig_genes), "'", collapse = ", ", sep = ""))
  
  
  # Summary of genes removed by neg corr
  message('\nThe following mutated genes were discarded due to negative correlations in expression and mutation profiles: ')
  
  message(paste("'", rownames(neg_cor_genes), "'", collapse = ", ", sep = ""))
  
  # Summary of genes removed by localization threshold
  excluded_genes <- genes_no_cor[!(rownames(genes_no_cor) %in% rownames(sig_genes)), ]
  
  message('\nThe following mutated genes were discarded due to a median q value below the input threshold: ')
  
  message(paste("'", rownames(excluded_genes), "'", collapse = ", ", sep = ""))

  
  # # Summary of genes removed by q value 
  # message('\nThe following mutated genes were discarded due to a nonsignificant q value (False Discovery Rate via BH procedure): ')
  # 
  # for (j in 1:num_complexes){
  #   message(paste('\n In complex', j, 'of', num_complexes, ':'))
  #   
  #   genes <- excluded_localization[[j]]
  #   
  #   if (length(genes) == 0){
  #     genes <-  'None'
  #   } 
  #   
  #   message(paste("'", genes, "'", collapse = ", ", sep = ""))
  # }
  # 
  # 
  # # Summary of genes removed because negative correlations
  # message('\nThe following mutated genes were discarded due to negative correlations between expression and mutation data: ')
  # 
  # for (k in 1:num_complexes){
  #   message(paste('\n In complex', k, 'of', num_complexes, ':'))
  #   
  #   genes <- excluded_neg_cor[[k]]
  #   
  #   if (length(genes) == 0){
  #     genes <-  'None'
  #   } 
  # 
  #   message(paste("'", genes, "'", collapse = ", ", sep = ""))
  # }
  # 
  
  TDAmut_object@significant_genes <- sig_genes
  
  return(TDAmut_object)
}
