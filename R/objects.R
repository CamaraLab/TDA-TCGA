#' Class to handle all information in the TDAmut pipeline
#'
#' @slot expression_table Expression data data.frame or matrix (create_TDAmut_object)
#' @slot mutation_table Mutation data table (create_TDAmut_object)
#' @slot nerve_complexes Topological representations as nerve complexes (compute_complexes)
#' @slot mapper_intervals Range of intervals used in TDAmapper (compute_complexes)
#' @slot mapper_percents Range of percent overlaps used in TDAmapper (compute_complexes)
#' @slot filter_embedding Embedding used to create topological representations in TDAmapper (compute_complexes)
#' @slot mutation_matrix Reformatted mutation data (create_TDAmut_object)
#' @slot nonsyn_mutations Reformatted nonsynonymous mutation data (create_TDAmut_object)
#' @slot syn_mutations Reformatted synonymous mutation data (create_TDAmut_object)
#' @slot mutational_load Total mutations per sample (create_TDAmut_object)
#' @slot mutational_load_localization p and q values indicating significance of mutational load localization across the parameter space (compute_complexes)
#' @slot correlations p and q values indicating significance of negative correlation between expression and mutation data of a given gene across the parameter space (filter_genes)
#' @slot filtered_genes Genes passing user-defined thresholds on mutational frequency, ratio of nonsynonymous / total mutations, and correlation q values (filter_genes)
#' @slot gene_scores p and q values indicating significance of localization of nonsynonymously mutated genes across the parameter space (compute_gene_localization)
#' @slot summary_matrix Matrix summarizing the results of the pipelinesin the computation (compute_gene_localization)
#'
#' @exportClass TDAmut
#'

TDAmut_object <- setClass(
  Class = 'TDAmut',
  slots = c(
    expression_table = 'ANY',
    mutation_table = 'data.frame',
    nerve_complexes = 'list',
    mapper_intervals = 'numeric',
    mapper_percents = 'numeric',
    filter_embedding = 'data.frame',
    mutation_matrix = 'data.frame',
    nonsyn_mutations = 'data.frame',
    syn_mutations = 'data.frame',
    mutational_load = 'numeric',
    mutational_load_localization = 'list',
    correlations = 'list',
    correlation_rhos = 'list',
    filtered_genes = 'character',
    gene_scores = 'list',
    summary_matrix = 'data.frame'
  )
)
