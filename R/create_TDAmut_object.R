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
#'
#' @export

# library(igraph)
# library(Matrix)
# library(ggplot2)
# library(dplyr)
# library(purrr)
# library(bioDist)
# library(dimRed)
# library(TDAmapper)
# library(umap)
# library(dbscan)
# library(RayleighSelection)
# library(cccd)
# library(maftools)

create_TDAmut_object <- function(exp_table, mut_table) {

  ######## INPUT AND CLEAN ########

  exp_table <- read.csv("/home/rstudio/documents/Messy_test_data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F)
  rownames(exp_table) <- substr(rownames(exp_table),1,16)
  exp_table <- exp_table[!(rownames(exp_table) %in% no_mut_data),]
  mut_table <- read.csv('/home/rstudio/documents/TDA-TCGA/Test_Data/LGG_Muts.txt',row.names=1,header=T,stringsAsFactors = F)
  # exp_table <- (read.csv(exp_table, row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))
  # mut_table <- read.csv(mut_table, row.names=1, header=T, stringsAsFactors = F, na.strings=c("NA","NaN", " ", "?")
  
  if(duplicated(rownames(exp_table))) {
    warning('Cleaning duplicated samples detected in expression data')
    exp_table <- exp_table[!duplicated(rownames(exp_table)),]
  }

  if(duplicated(colnames(exp_table))) {
    warning('Cleaning duplicated genes detected in expression data')
    exp_table <- exp_table[,!duplicated(colnames(exp_table))]
  }

  if(all(unique(mut_table$Sample)) %in% rownames(exp_table)) {
    message("Matching samples between mutation and expression data")
  } else {
    no_mut_data <- !(rownames(exp_table) %in% unique(mut_table$Sample))
    warning(paste0("Removing expression samples with no mutation data: ", no_mut_data))
    exp_table <- exp_table[!(rownames(exp_table) %in% no_mut_data),]
  }

  TDAmut_object <- new(
    Class = 'TDAmut',
    expression_table = exp_table
    mutation_table = mut_table
  )
  
  return(TDAmut_object)
}