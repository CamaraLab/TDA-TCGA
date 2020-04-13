#' Class to handle all information in TDA.TCGA pipeline
#' 
#' @slot nerve_complexes
#' @slot mapper_intervals
#' @slot mapper_percents
#' @slot expression_table
#' @slot mutation_table
#' @slot filter_embedding
#' @slot nonsyn_mutations
#' @slot syn_mutations
#' @slot mutational_load
#' @slot mutational_load_localization
#' @slot gene_localization
#' @slot significant_genes
#' 
#' @exportClass TDAmut
#' 
#' 

TDAmut_object <- setClass(
  Class = 'TDAmut',
  slots = c(
    nerve_complexes = 'list',
    mapper_intervals = 'numeric',
    mapper_percents = 'numeric',
    expression_table = 'ANY', 
    mutation_table = 'data.frame',
    filter_embedding = 'data.frame',
    nonsyn_mutations = 'data.frame',
    syn_mutations = 'data.frame',
    mutational_load = 'numeric',
    mutational_load_localization = 'list',
    gene_localization = 'list',
    significant_genes = 'data.frame'
  )
)