#' Computes localization of nonsynonymously mutated genes passing thresholds of filter_genes
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes
#' @param neg_correlations_threshold q value theshold to remove genes displaying correlations between mutation and expression data. Genes with a median q value across all complexes exceeding this threshold are not considered. If NULL, negative correlations will be ignored. By default is 0.8.
#' @param num_permutations number of permutations used to create null distribution when assessing localization of mutations in nerve complexes. By default is 10,000
#' @param num_cores number of cores to be used in computation. By default is 1.
#' @param seed integer specifying the seed used to initialize the generator of permutations. By default is 121
#'
#' @return Returns TDAmut object populated with Laplacian scores, p values, and q values for localization of genes passing thresholds of filter_genes
#'
#' @export

compute_gene_localization <- function(TDAmut_object, negative_correlations_threshold = 0.8, num_permutations = 10000, num_cores = 5, seed = 121){

  if (is_empty(TDAmut_object@nerve_complexes) || is_empty(TDAmut_object@nonsyn_mutations) || is_empty(TDAmut_object@filtered_genes)){
    stop('Run compute_complexes, compute_mut_load, and filter_genes first to populate object with nerve complexes, formatted mutation data, and candidate genes')
  }

  nerve_complexes <- TDAmut_object@nerve_complexes
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  genes <- TDAmut_object@filtered_genes

  if (!is.null(negative_correlations_threshold)){

    if(is_empty(TDAmut_object@negative_correlations)){
      stop('Run filter_genes with negative_correlations = TRUE first')
    }

    neg_cor_list <- TDAmut_object@negative_correlations

    # Rearranging scores by gene instead of by nerve complex
    collapsed <- bind_rows(neg_cor_list)
    neg_cor_q_bygene <- matrix(collapsed$q, length(genes), length(nerve_complexes), dimnames = list(unique(collapsed$Gene), NULL), byrow = FALSE) %>% as.data.frame # genes x complexes
    neg_cor_genes <- neg_cor_q_bygene[apply(neg_cor_q_bygene, 1, function(x) any(median(x) > negative_correlations_threshold)), ]

    genes <- genes[!(genes %in% rownames(neg_cor_genes))] # Not considering genes with neg correlation q val above threshold

    message('Not considering the following ', nrow(neg_cor_genes), ' genes displaying negative correlations between mutation and expression data: ', paste("'", rownames(neg_cor_genes), "'", collapse = ", ", sep = ""))

  }

  ######## ASSESSING LOCALIZATION OF GENES IN NERVE COMPLEXES ########

  nonsyn_features <- t(nonsyn_mat[ , genes]) %>% as.data.frame
  num_complexes <- length(nerve_complexes)
  gene_scores <- vector('list', num_complexes)

  for(i in 1:num_complexes){
    message(paste("Assessing genes in complex", i, "of", num_complexes))
    gene_scores[[i]] <- rayleigh_selection(nerve_complexes[[i]], nonsyn_features, num_perm = num_permutations, num_cores = num_cores, one_forms = FALSE, seed = seed)
    gene_scores[[i]] <- data.frame("Gene"= rownames(nonsyn_features), "R" = sapply(gene_scores[i], '[', 1), "p" = sapply(gene_scores[i], '[', 2), "q" = sapply(gene_scores[i], '[', 3))
  }

  TDAmut_object@gene_scores <- gene_scores
  return(TDAmut_object)

}
