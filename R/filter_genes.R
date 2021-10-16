#' Filter genes considered for compute_localization by mutation frequency, fraction of nonsynonymous mutations, and negative correlations between expression and mutation data
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes
#' @param freq_threshold threshold mutation frequency. Genes below this value are not considered in analysis. By default is 0.02
#' @param top_nonsyn_fraction number of genes to keep with greatest nonsyn/nonsyn+syn fraction. By default is 350.
#' @param upper_correlations_threshold upper q value theshold to keep genes displaying positive correlations between mutation and expression data. Genes with a median q value across all complexes exceeding this threshold are considered for further analysis. By default is 0.9.
#' @param lower_correlations_threshold lower q value theshold to keep genes displaying negative correlations between mutation and expression data. Genes with a median q value across all complexes below this threshold are considered for further analysis. By default is 1e-4.
#'
#' @return Returns a TDAmut object populated with filtered genes to consider in compute_gene_localization. Optionally returns p and q values quantifying negative correlations between expression and mutation profiles of filtered genes.
#'
#' @export

filter_genes <- function(TDAmut_object, freq_threshold = 0.02, top_nonsyn_fraction = 350, upper_correlations_threshold = 0.9, lower_correlations_threshold = 1e-4) {

  if (is_empty(TDAmut_object@nerve_complexes)){
    stop('Run compute_complexes first to populate object with nerve complexes')
  }

  nerve_complexes <- TDAmut_object@nerve_complexes
  exp_table <- TDAmut_object@expression_table
  mut_mat <- TDAmut_object@mutation_matrix
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  syn_mat <- TDAmut_object@syn_mutations

  ######## APPLYING FREQUENCY AND NONSYNONYMOUS FRACTION THRESHOLD ########

  mut_mat <- mut_mat[order(rownames(mut_mat)), order(colnames(mut_mat))]
  mut_bin <- ifelse(mut_mat > 0, 1, 0) %>% as.data.frame

  nonsyn_bin <- ifelse(nonsyn_mat > 0, 1, 0) %>% as.data.frame
  freq <- colMeans(nonsyn_bin) %>% as.data.frame
  nonsyn_mat <- nonsyn_mat[ , freq > freq_threshold]

  if (ncol(nonsyn_mat) >= top_nonsyn_fraction){
    nonsyn_frac <- colSums(nonsyn_mat) / colSums(nonsyn_mat + syn_mat[ , colnames(nonsyn_mat)])
    nonsyn_mat <- nonsyn_mat[ , order(nonsyn_frac, decreasing = TRUE)][ , 1:top_nonsyn_fraction]
  }

  genes <- colnames(nonsyn_mat)

  if (!all(genes %in% colnames(exp_table))){
    missing_genes <- genes[!genes %in% colnames(exp_table)]
    genes <- genes[!genes %in% missing_genes]
    message('Not considering the following genes with no expression data: ', paste("'", missing_genes, "'", collapse = ", ", sep = ""))
  }

  ######## COMPUTING CORRELATIONS ########

  num_complexes <- length(nerve_complexes)
  cor_list <- replicate(num_complexes, data.frame('Gene' = NULL, 'p' = NULL), simplify = F) %>% set_names(1:num_complexes)
  rho_list <- replicate(num_complexes, data.frame('Gene' = NULL, 'rho' = NULL), simplify = F) %>% set_names(1:num_complexes)

  # auxiliary function, used below for computing the mean value of some feature over the network
  push <- function(g2, lo, pushforward) {
    return(as.numeric(lapply(g2$points_in_vertex, function(idx) pushforward(lo[idx]))))
  }

  count = 0
  for (complex in nerve_complexes) {
    count = count + 1

    message(paste("Computing Correlations Between Gene Expression and Gene Mutation In Complex", count, "of", num_complexes))

    ps <- rep(NA, length(genes)) %>% as.list
    rhos <- rep(NA, length(genes)) %>% as.list

    gene_count = 1
    for (gene in genes) {
      mut_mean <- as.numeric(push(g2=complex, lo=mut_mat[, colnames(mut_mat) == gene], pushforward = mean))
      exp_mean <- as.numeric(push(g2=complex, lo=exp_table[, colnames(exp_table) == gene], pushforward = mean))

      spearman_cor <- tryCatch(
        {
          suppressWarnings(cor.test(mut_mean, exp_mean, method="spearman", alternative="less")) # alternative = “less”, so lower p-val means more negatively correlated
        },
        error=function(cond) {
          message(paste("Error in complex", count, "with gene", gene))
          message(cond)
          message("\n")
          return(1)
        })

      ps[[gene_count]] <- spearman_cor$p.value
      rhos[[gene_count]] <- spearman_cor$estimate
      gene_count = gene_count + 1
    }

    cor_list[[count]] <- data.frame('Gene' = genes, 'p' = unlist(ps))
    rho_list[[count]] <- data.frame('Gene' = genes, 'rho' = unlist(rhos))
  }

  cor_list <- lapply(cor_list, function(x) {
    q <- p.adjust(x$p, method = 'BH')
    x <- cbind(x, q)
  })

  TDAmut_object@correlations <- cor_list
  TDAmut_object@correlation_rhos <- rho_list

  ######## FILTERING OUT GENES BY CORRELATION ########

  if (!(is.null(upper_correlations_threshold)| is.null(lower_correlations_threshold))){
    # Rearranging scores by gene instead of by nerve complex
    collapsed <- bind_rows(cor_list)
    cor_q_bygene <- matrix(collapsed$q, length(genes), length(nerve_complexes), dimnames = list(unique(collapsed$Gene), NULL), byrow = FALSE) %>% as.data.frame # genes x complexes
    weak_cor_genes <- cor_q_bygene[apply(cor_q_bygene, 1, function(x) (median(x) < upper_correlations_threshold & median(x) > lower_correlations_threshold)), ]

    genes <- genes[!(genes %in% rownames(weak_cor_genes))] # Not considering genes with correlation q val in between upper threshold and lower threshold

    message('Not considering the following ', nrow(weak_cor_genes), ' genes not displaying strong correlations/anticorrelations between mutation and expression data: ', paste("'", rownames(weak_cor_genes), "'", collapse = ", ", sep = ""))
  }

  TDAmut_object@filtered_genes <- genes

  return(TDAmut_object)

}

