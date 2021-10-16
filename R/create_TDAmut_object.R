#' Creates object of class TDAmut to use throughout pipeline
#'
#' @param exp_table transcriptomic data from cohort with format samples x genes
#' @param mut_table mutation data from cohort as a data.frame with columns Sample, Gene, Mutation, Type
#' @param samples_to_use optional list of samples names to excluded from the analysis
#'
#' @import igraph
#' @import Matrix
#' @import ggplot2
#' @import dplyr
#' @import purrr
#' @import bioDist
#' @import dimRed
#' @import TDAmapper
#' @import umap
#' @import dbscan
#' @import RayleighSelection
#' @import cccd
#' @import fields
#' @import reshape2
#' @import Rcpp
#' @import parallel
#'
#' @return Creates a TDAmut object populated with expression data, mutation data, and data frames of nonsynonymous mutations, synonymous mutations, and mutational load of samples
#'
#' @export

create_TDAmut_object <- function(exp_table, mut_table, samples_to_exclude = NULL) {

  ######## INPUT AND CLEAN ########

  # exp_table <- read.csv("/home/rstudio/documents/TDA-TCGA/data/LGG_Full_TPM.csv", row.names = 1, header = T, stringsAsFactors = F)
  # rownames(exp_table) <- substr(rownames(exp_table), 1, 16)
  # mut_table <- read.csv('/home/rstudio/documents/TDA-TCGA/data/LGG_Muts.txt', row.names = 1, header = T, stringsAsFactors = F)
  exp_table <- read.csv(exp_table, row.names = 1, header = T, stringsAsFactors = F, na.strings=c("NA","NaN", " ", "?"))
  mut_table <- read.csv(mut_table, row.names = 1, header = T, stringsAsFactors = F, na.strings=c("NA","NaN", " ", "?"))

  # remove genes for which expression is 0 (sum(vector)==0)
  if(any(colSums((exp_table)) == 0)) {
    num_genes_removed <- length(colnames(exp_table)[colSums(exp_table) == 0])
    exp_table <- exp_table[, colSums(exp_table) != 0]
    message('Removed ', paste(num_genes_removed), ' genes with no expression data')
  }

  # remove rows from the mut table which don't have a corresponding gene
  if(any(mut_table$Gene == '.')) {
    num_rows_removed <- length(rownames(mut_table)[mut_table$Gene == '.'])
    mut_table <- mut_table[mut_table$Gene != '.', ]
    message('Removed ', paste(num_rows_removed), ' mutation entries with no corresponding gene')
  }

  if(any(duplicated(rownames(exp_table)))) {
    exp_table <- exp_table[!duplicated(rownames(exp_table)), ]
    message('Removed duplicated samples detected in expression data')
  }

  if(any(duplicated(colnames(exp_table)))) {
    exp_table <- exp_table[ , !duplicated(colnames(exp_table))]
    message('Removed duplicated genes detected in expression data')
  }

  if(!(all(unique(mut_table$Sample) %in% rownames(exp_table)))) {
    no_exp_data <- mut_table$Sample[!(unique(mut_table$Sample) %in% rownames(exp_table))]
    mut_table <- mut_table[!(mut_table$Sample %in% no_exp_data), ]
    message('Removed samples in mutation data not in expression data: ', paste("'",no_exp_data,"'",collapse=", ",sep=""))
  }
  else if(!(all(rownames(exp_table) %in% unique(mut_table$Sample)))) {
    no_mut_data <- rownames(exp_table[!(rownames(exp_table) %in% unique(mut_table$Sample)), ])
    exp_table <- exp_table[!(rownames(exp_table) %in% no_mut_data), ]
    message('Removed samples in expression data not in mutation data: ', paste("'",no_mut_data,"'",collapse=", ",sep=""))
  }
  else {
    message("Samples match between expression and mutation data")
  }

  if (!all(unique(mut_table$Gene) %in% colnames(exp_table))){
    missing_genes_exp <- unique(mut_table$Gene[!(mut_table$Gene %in% colnames(exp_table))])
    message('The following genes have mutation data but no expression data. They will not be considered for filtering of negative correlations later in the TDAmut pipeline: ',
            paste("'", missing_genes_exp, "'", collapse = ", ", sep = ""))
  }

  ######## OPTIONAL SUBSETTING OF SAMPLES ########

  if (!is.null(samples_to_exclude)) {
    exp_table = exp_table[!(rownames(exp_table) %in% samples_to_exclude),]
    mut_table = mut_table[!(mut_table$Sample %in% samples_to_exclude),]
  }

  ######## CONSOLIDATE MUTATION TABLE INTO DATA FRAMES ########

  split_mut_data <- function(TDAmut_object) {

    samples <- sort(unique(mut_table$Sample))
    gene_names <- sort(unique(mut_table$Gene))

    # Reformatting table into matrix
    mut_mat <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
    dimnames(mut_mat) <- list(samples,gene_names)
    t_mut <- with(mut_table,table(Sample,Gene))
    mut_mat[rownames(t_mut),colnames(t_mut)] <- t_mut
    mut_mat <- mut_mat[order(rownames(mut_mat)), order(colnames(mut_mat))]

    # Assigning syn and nonsyn mutations
    if(is.null(nonsyn_muts) || is.null(syn_muts)){
      nonsyn_type <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
                       "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins",
                       "In_Frame_Ins", "Nonstop_Mutation")

      nonsyn_muts <- mut_table[mut_table$Type %in% nonsyn_type,]
      syn_muts <- mut_table[!(mut_table$Type %in% nonsyn_type),]
    } else {
      nonsyn_muts <- read.csv(nonsyn_muts, row.names = 1, header = T, stringsAsFactors = F, na.strings=c("NA","NaN", " ", "?"))
      syn_muts <- read.csv(syn_muts, row.names = 1, header = T, stringsAsFactors = F, na.strings=c("NA","NaN", " ", "?"))
    }

    # Reformatting nonsynonymous and synonymous mutations into matrices
    nonsyn_mat <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
    dimnames(nonsyn_mat) <- list(samples,gene_names)
    syn_mat <- nonsyn_mat

    t_nonsyn <- with(nonsyn_muts,table(Sample,Gene))
    nonsyn_mat[rownames(t_nonsyn),colnames(t_nonsyn)] <- t_nonsyn
    nonsyn_mat <- nonsyn_mat[order(rownames(nonsyn_mat)), order(colnames(nonsyn_mat))]
    t_syn <- with(syn_muts,table(Sample,Gene))
    syn_mat[rownames(t_syn),colnames(t_syn)] <- t_syn
    syn_mat <- syn_mat[order(rownames(syn_mat)), order(colnames(syn_mat))]

    # re-order the three matrices so they have the same ordering as the exp matrix
    mut_mat <- mut_mat[rownames(exp_table), ]
    nonsyn_mat <- nonsyn_mat[rownames(exp_table), ]
    syn_mat <- syn_mat[rownames(exp_table), ]

    # removing NA rows (these exist because of samples that get removed)
    mut_mat <- mut_mat[!is.na(nonsyn_mat[, 1]), ]
    syn_mat <- syn_mat[!is.na(nonsyn_mat[, 1]), ]
    exp_table <- exp_table[!is.na(nonsyn_mat[, 1]), ]
    nonsyn_mat <- nonsyn_mat[!is.na(nonsyn_mat[, 1]), ]

    # removing duplicated samples (removes any row associated with a sample of the form/name <sample>.<number>; this only occurs if there already exists a row/sample name of the form <sample> in the table)
    # we need the function because if the below call to grep returns an empty vector, [-<vector>, ] will just delete the entire matrix :(
    remove_duplicated_samples <- function(matrix, rows_to_remove) {
      if (length(rows_to_remove) > 0) {
        return(matrix[-rows_to_remove, ])
      } else {
        return(matrix)
      }
    }
    nonsyn_mat <- remove_duplicated_samples(nonsyn_mat, grep(".\\..", rownames(nonsyn_mat)))
    syn_mat <- remove_duplicated_samples(syn_mat, grep(".\\..", rownames(syn_mat)))
    mut_mat <- remove_duplicated_samples(mut_mat, grep(".\\..", rownames(mut_mat)))
    exp_table <- remove_duplicated_samples(exp_table, grep(".\\..", rownames(exp_table)))

    TDAmut_object@mutation_matrix <- mut_mat
    TDAmut_object@nonsyn_mutations <- nonsyn_mat
    TDAmut_object@syn_mutations <- syn_mat
    TDAmut_object@expression_table <- exp_table

    TDAmut_object@mutational_load <- rowSums(syn_mat + nonsyn_mat)

    return(TDAmut_object)
  }

  ######## CREATING AND POPULATING OBJECT ########

  TDAmut_object <- new(
    Class = 'TDAmut',
    expression_table = exp_table,
    mutation_table = mut_table
  )

  TDAmut_object <- split_mut_data(TDAmut_object)

  return(TDAmut_object)
}
