#' Class to handle all information in TDA.TCGA pipeline
#' 
#' @slot nerve_complexes
#' @slot mapper_intervals
#' @slot mapper_percents
#' @slot expression_table
#' @slot mutation_table
#' @slot mutation_matrix
#' @slot filter_embedding
#' @slot nonsyn_mutations
#' @slot syn_mutations
#' @slot mutational_load
#' @slot min_mutated_samples
#' @slot mutational_load_localization
#' @slot gene_scores
#' @slot significant_genes
#' @slot negative_correlations
#' @slot filtered_genes
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
    min_mutated_samples = 'ANY',
    mutational_load_localization = 'list',
    gene_scores = 'list',
    significant_genes = 'list',
    negative_correlations = 'list',
    filtered_genes = 'character'
  )
)
