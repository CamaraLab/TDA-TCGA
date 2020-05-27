#' Class to handle all information in the TDAmut pipeline
#'
#' @slot expression_table Expression data data.frame or matrix (create_TDAmut_object)
#' @slot mutation_table Mutation data table (create_TDAmut_object)
#' @slot nerve_complexes Topological representations as nerve complexes (compute_complexes)
#' @slot mapper_intervals Range of intervals used in TDAmapper (compute_complexes)
#' @slot mapper_percents Range of percent overlaps used in TDAmapper (compute_complexes)
#' @slot filter_embedding Embedding used to create topological representations in TDAmapper (compute_complexes)
#' @slot mutation_matrix Reformatted mutation data (compute_mut_load)
#' @slot nonsyn_mutations Reformatted nonsynonymous mutation data (compute_mut_load)
#' @slot syn_mutations Reformatted synonymous mutation data (compute_mut_load)
#' @slot mutational_load Total mutations per sample (compute_mut_load)
#' @slot min_mutated_samples Samples with mutational load below a user-defined threshold (compute_mut_load)
#' @slot mutational_load_localization p and q values indicating significance of mutational load localization across the parameter space (compute_mut_load)
#' @slot negative_correlations p and q values indicating significance of correlation between expression and mutation data of a given gene across the parameter space (filter_genes)
#' @slot filtered_genes Genes passing user-defined thresholds on mutational frequency and ratio of nonsynonymous / total mutations
#' @slot gene_scores p and q values indicating significance of localization of nonsynonymously mutated genes across the parameter space (compute_gene_localization)
#' @slot significant_genes Genes passing user-defined median q value threshold across the parameter space (identify_significant_mutations)
#'
#' @exportClass TDAmut
#'

TDAmut_object <- setClass(
  Class = 'TDAmut',
  slots = c(
    nerve_complexes = 'list',
    mapper_intervals = 'numeric',
    mapper_percents = 'numeric',
    expression_table = 'ANY',
    mutation_table = 'data.frame',
    mutation_matrix = 'data.frame',
    filter_embedding = 'data.frame',
    nonsyn_mutations = 'data.frame',
    syn_mutations = 'data.frame',
    mutational_load = 'numeric',
    min_mutated_samples = 'character',
    mutational_load_localization = 'list',
    gene_scores = 'list',
    significant_genes = 'list',
    negative_correlations = 'list',
    filtered_genes = 'character'
  )
)
