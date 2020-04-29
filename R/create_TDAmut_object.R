#' Creates object of class TDAmut to use throughout pipeline
#'
#' @param exp_table transcriptomic data from cohort with format samples x genes
#' @param mut_table mutation data from cohort as a data.frame with columns Sample, Gene, Mutation, Type
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
#' @import maftools
#' @import fields
#' @import reshape2
#'
#' @export

library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)
library(purrr)
library(bioDist)
library(dimRed)
library(TDAmapper)
library(umap)
library(dbscan)
library(RayleighSelection)
library(cccd)
library(maftools)
library(fields)
library(reshape2)

create_TDAmut_object <- function(exp_table, mut_table) {

  ######## INPUT AND CLEAN ########

  exp_table <- read.csv("/home/rstudio/documents/TDA-TCGA/data/LGG_Full_TPM_matrix.csv", row.names = 1, header = T, stringsAsFactors = F)
  rownames(exp_table) <- substr(rownames(exp_table), 1, 16)
  #exp_table <- exp_table[!(rownames(exp_table) %in% no_mut_data),]
  mut_table <- read.csv('/home/rstudio/documents/TDA-TCGA/data/LGG_Muts.txt', row.names = 1, header = T, stringsAsFactors = F)
  # exp_table <- (read.csv(exp_table, row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))
  # mut_table <- read.csv(mut_table, row.names=1, header=T, stringsAsFactors = F, na.strings=c("NA","NaN", " ", "?")
  
  
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
    missing_genes_exp <- unique(mut_table$Gene[unique(mut_table$Gene) %in% colnames(exp_table)])
    warning('The following genes have mutation data but no expression data, 
            which limits the optional filtering of negative correlations later in the TDAmut pipeline: ',
            paste("'", missing_genes_exp, "'", collapse = ", ", sep = ""))
  }

  ######## CREATING AND POPULATING OBJECT ########
  
  TDAmut_object <- new(
    Class = 'TDAmut',
    expression_table = exp_table,
    mutation_table = mut_table
  )
  
  return(TDAmut_object)
}