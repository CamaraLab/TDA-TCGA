#' Controls for negative correlations between mutations rates and expression of a gene
#' 
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, and nerve complexes
#' @param significant_genes if TRUE, only considers genes output by identify_significant_mutations. By default is FALSE
#' @param num_permutations number of permutations used to create null distribution. By default is 2000
#' 
#' @return Returns TDAmut object populated with p and q values indicating correlations between expression and mutation data of genes
#' 
#' @import philentropy
#' 
#' @export
#' 
filter_negative_correlations <- function(TDAmut_object, significant_genes = FALSE, num_permutations = 2000){
  
  if (significant_genes == TRUE){
    gene_scores_table <- TDAmut_object@significant_genes
    genes <- unique(unlist(sapply(gene_scores,'[[', 1))) %>% as.character
  }
  else{
    gene_scores_table <- TDAmut_object@gene_scores 
    genes <- (gene_scores_table[[1]])$Gene %>% as.character
  }
  
  if (is.null(TDAmut_object@nerve_complexes) || is.null(TDAmut_object@mutation_matrix) || is.null(genes)){
    stop('Please run compute_complexes, compute_mut_load, and assess_mutations first to populate object with nerve complexes, formatted mutation data, and gene localization scores')
  }
  
  nerve_complexes <- TDAmut_object@nerve_complexes
  exp_table <- TDAmut_object@expression_table
  mut_mat <- TDAmut_object@mutation_matrix
  
  mut_mat <- mut_mat[order(rownames(mut_mat)), order(colnames(mut_mat))]
  mut_bin <- ifelse(mut_mat > 0, 1, 0) %>% as.data.frame
  
  if (!all(genes %in% colnames(exp_table))){
    missing_genes <- genes[!genes %in% colnames(exp_table)]
    genes <- genes[!genes %in% missing_genes]
    message('Not considering the following genes with no expression data: ', paste("'", missing_genes, "'", collapse = ", ", sep = ""))
  }
  
  vertices_in_graphs <- lapply(nerve_complexes, function(x) {
    x <- x$points_in_vertex
    return(x)
  })

  graph_list <- replicate(length(nerve_complexes), data.frame('Gene' = NULL, 'p' = NULL), simplify = F)
  count = 0
  
  for (graph in vertices_in_graphs){ # iterates through grid of complexes
    count = count + 1
    message(paste('Evaluating genes in complex', count, 'of', length(nerve_complexes)))
    for (gene in genes){ # iterates through genes of interest
      
      mut_freq <- sapply(graph, function(x){ # ei
        mean(mut_bin[x, gene])
        })
      
      avg_exp <- sapply(graph, function(x){ # ri
        mean(exp_table[x, gene])
      })
      
      mut_perms <- sapply(1:num_permutations, function(x){ # returns table with rows = vertices, columns = permutations
        perm <- mut_bin[sample(nrow(mut_bin)), ]
        sapply(graph, function(y) mean(perm[y, gene]))
      })
      
      mut_freq_norm <- mut_freq / sum(mut_freq)
      mut_perms_norm <- sweep(mut_perms, 2, colSums(mut_perms), `/`)
      avg_exp_norm <- avg_exp / sum(avg_exp)
      
      jsdiv <- suppressMessages(JSD(rbind(avg_exp_norm, mut_freq_norm)))
      jsdiv_perm <- suppressMessages(JSD(rbind(avg_exp_norm, t(mut_perms_norm))))
      
      p <- sum(tail(jsdiv_perm[ , 1], -1) < jsdiv) / num_permutations
      
      graph_list[[count]] <- rbind(graph_list[[count]], data.frame(gene, p))
    }
    
    q <- p.adjust(graph_list[[count]]$p, method = 'BH')
    graph_list[[count]] <- cbind(graph_list[[count]], q)
    
  }

  return(graph_list)
  # TDAmut_object@negative_correlations <- graph_list
  # return(TDAmut_object)
  
  

  
}



