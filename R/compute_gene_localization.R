#' Computes localization of nonsynonymously mutated genes passing thresholds of filter_genes
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes
#' @param num_permutations number of permutations used to create null distribution when assessing localization of mutations in nerve complexes. By default is 10,000
#' @param num_cores number of cores to be used in computation. By default is 1.
#' @param seed integer specifying the seed used to initialize the generator of permutations. By default is 121
#'
#' @return Returns TDAmut object populated with Laplacian scores, p values, and q values for localization of genes passing thresholds of filter_genes. Outputs genes which are significantly localized.
#'
#' @export

compute_gene_localization <- function(TDAmut_object, num_permutations = 10000, num_cores = 1, seed = 121){

  if (is_empty(TDAmut_object@filtered_genes)){
    stop('Run filter_genes first to populate object with candidate genes')
  }

  ######## ASSESSING LOCALIZATION OF GENES IN NERVE COMPLEXES ########

  nerve_complexes <- TDAmut_object@nerve_complexes
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  genes <- TDAmut_object@filtered_genes

  nonsyn_features <- t(nonsyn_mat[ , genes]) %>% as.data.frame
  num_complexes <- length(nerve_complexes)
  gene_scores <- vector('list', num_complexes)

  # mutational load is a covariate
  mutload <- TDAmut_object@mutational_load
  mutload_mat <- t(mutload) %>% as.data.frame
  message(paste("Computing Gene Scores"))

  # supressing the "Inf values in sample" warning, since it does not indicate problems in the computation. it also does nt appear often
  suppressWarnings(gene_scores <- rayleigh_selection(nerve_complexes, nonsyn_features, num_perms = num_permutations, num_cores = num_cores, one_forms = FALSE, seed = seed, covariates = mutload_mat, optimize.p = "perm"))

  # the raw output from rayleigh_selection
  TDAmut_object@gene_scores <- gene_scores

  ######## REFORMATING LOCALIZATION OUTPUT #######

  message(paste("Gene Scores Computed. Formatting Output"))

  summary_mat <- gene_scores[order(gene_scores$combined.p0), colnames(gene_scores)=="combined.p0" | colnames(gene_scores)=="q0"]

  # compute sumamry statistics on expression for each gene
  exp_table <- TDAmut_object@expression_table
  first_quartile_exp <- lapply(rownames(summary_mat), function(x) unname(quantile(exp_table[, colnames(exp_table) == x], 0.25)))
  median_exp <- lapply(rownames(summary_mat), function(x) median(exp_table[, colnames(exp_table) == x]))
  third_quartile_exp <- lapply(rownames(summary_mat), function(x) unname(quantile(exp_table[, colnames(exp_table) == x], 0.75)))

  summary_mat$expression_first_quartile <- unlist(first_quartile_exp)
  summary_mat$expression_median <- unlist(median_exp)
  summary_mat$expression_third_quartile <- unlist(third_quartile_exp)

  # extract summary statistics related to the negative correlation computations for each gene
  cor_list <- TDAmut_object@correlations
  collapsed <- bind_rows(cor_list)
  collapsed <- collapsed[collapsed$Gene %in% genes,]

  cor_p_bygene <- matrix(collapsed$p, length(genes), length(nerve_complexes), dimnames = list(unique(collapsed$Gene), NULL), byrow = FALSE) %>% as.data.frame # genes x complexes
  cor_median_ps <- apply(cor_p_bygene, 1, median)
  cor_median_ps <- cor_median_ps[rownames(summary_mat)]

  cor_q_bygene <- matrix(collapsed$q, length(genes), length(nerve_complexes), dimnames = list(unique(collapsed$Gene), NULL), byrow = FALSE) %>% as.data.frame # genes x complexes
  cor_median_qs <- apply(cor_q_bygene, 1, median)
  cor_median_qs <- cor_median_qs[rownames(summary_mat)]

  rhos <- TDAmut_object@correlation_rhos
  collapsed_rhos <- bind_rows(rhos)
  collapsed_rhos <- collapsed_rhos[collapsed_rhos$Gene %in% genes,]

  cor_rho_bygene <- matrix(collapsed_rhos$rho, length(genes), length(nerve_complexes), dimnames = list(unique(collapsed_rhos$Gene), NULL), byrow = FALSE) %>% as.data.frame
  cor_median_rhos <- apply(cor_rho_bygene, 1, median)
  cor_median_rhos <- cor_median_rhos[rownames(summary_mat)]

  summary_mat$cor_median_p <- cor_median_ps
  summary_mat$cor_median_q <- cor_median_qs
  summary_mat$cor_median_rho <- cor_median_rhos

  # number of patients with a nonsynonmous mutation in each gene
  nonsyn_count <- lapply(rownames(summary_mat), function(x) sum(nonsyn_mat[, colnames(nonsyn_mat) == x] > 0))
  summary_mat$nonsyn_count <- unlist(nonsyn_count)

  TDAmut_object@summary_matrix <- summary_mat

  return(TDAmut_object)

}
