#' Compute each gene's Laplacian score, p value, and q value via BH procedure for controlling false discovery rate across grid of nerve complexes
#' 
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes
#' @param freq_threshold threshold mutation frequency. Genes below this value are not considered in analysis. By default is 0.02
#' @param top_nonsyn_fraction number of genes to keep with greatest nonsyn/nonsyn+syn fraction. By default is 350.
#' @param negative_correlations if TRUE, assesses negative correlations between mutations rates and expression of a gene. By default is TRUE.
#' @param num_permutations_negative_correlations if negative_correlations = TRUE, number of permutations used to create null distribution when identifying negative correlations between expression and mutation data. By default is 2000
#' @param num_permutations_gene_localization number of permutations used to create null distribution when assessing localization of mutations in nerve complexes. By default is 10,000
#' @param num_cores number of cores to be used in computation. By default is 1.
#' @param seed integer specifying the seed used to initialize the generator of permutations. By default is 121
#' 
#' @import philentropy
#' @import RayleighSelection
#' 
#' @return Returns TDAmut object populated with filtered genes' Laplacian scores, p values, and q values for localization and negative correlations
#' 
#' @export

filter_genes <- function(TDAmut_object, freq_threshold = 0.01, top_nonsyn_fraction = 350, negative_correlations = TRUE, num_permutations_negative_correlation = 2000, num_permutations_gene_localization = 10000, num_cores = 5, seed = 121){
  
  if (is_empty(TDAmut_object@nerve_complexes) || is_empty(TDAmut_object@nonsyn_mutations)){
    stop('Run compute_complexes and compute_mut_load first to populate object with nerve complexes and formatted mutation data')
  }
  
  nerve_complexes <- TDAmut_object@nerve_complexes
  exp_table <- TDAmut_object@expression_table
  mut_mat <- TDAmut_object@mutation_matrix
  nonsyn_mat <- TDAmut_object@nonsyn_mutations
  syn_mat <- TDAmut_object@syn_mutations
  
  mut_mat <- mut_mat[order(rownames(mut_mat)), order(colnames(mut_mat))]
  mut_bin <- ifelse(mut_mat > 0, 1, 0) %>% as.data.frame
  
  ######## APPLYING FREQUENCY THRESHOLD ########
  
  nonsyn_bin <- ifelse(nonsyn_mat > 0, 1, 0) %>% as.data.frame
  freq <- colMeans(nonsyn_bin) %>% as.data.frame
  nonsyn_mat <- nonsyn_mat[ , freq > freq_threshold]
  
  ######## APPLYING NONSYNONYMOUS FRACTION THRESHOLD ########
  
  if (ncol(nonsyn_mat) >= top_nonsyn_fraction){
    nonsyn_frac <- colSums(nonsyn_mat) / colSums(nonsyn_mat + syn_mat[ , colnames(nonsyn_mat)])
    nonsyn_mat <- nonsyn_mat[ , order(nonsyn_frac, decreasing = TRUE)][ , 1:top_nonsyn_fraction]
  }
  
  ######## IDENTIFYING NEGATIVE CORRELATIONS ########
  
  if (negative_correlations == TRUE){
    
    genes <- colnames(nonsyn_mat)
    if (!all(genes %in% colnames(exp_table))){
      missing_genes <- genes[!genes %in% colnames(exp_table)]
      genes <- genes[!genes %in% missing_genes]
      message('Not considering the following genes with no expression data: ', paste("'", missing_genes, "'", collapse = ", ", sep = ""))
    }
    
    genes <- LGG_drivers ##DELETE
    
    vertices_in_graphs <- lapply(nerve_complexes, function(x) {
      x <- x$points_in_vertex
      return(x)
    })
    
    neg_cor_list <- replicate(length(nerve_complexes), data.frame('Gene' = NULL, 'p' = NULL), simplify = F)
    count = 0
    RNGkind("L'Ecuyer-CMRG") # For reproducible random number generation with forking parallelization
    set.seed(seed)
    
    
    for (graph in vertices_in_graphs){ 
      count = count + 1
      message(paste('Evaluating negative correlations in complex', count, 'of', length(nerve_complexes)))
      for (gene in genes){
        
        message('CHECKING: ', gene)
        
        # mut_freq <- sapply(graph, function(x){ # ei
        #   mean(mut_bin[x, gene])
        # })
        # 
        # avg_exp <- sapply(graph, function(x){ # ri
        #   mean(exp_table[x, gene])
        # })
        # 
        # mut_perms <- sapply(1:num_permutations_negative_correlation, function(x){ # returns table with rows = vertices, columns = permutations
        #   perm <- mut_bin[sample(nrow(mut_bin)), ]
        #   sapply(graph, function(y) mean(perm[y, gene]))
        # })
        
        mut_freq <- unlist(mclapply(graph, mc.cores = num_cores, function(x){ # ei
          mean(mut_bin[x, gene])
        }))
        
        avg_exp <- unlist(mclapply(graph, mc.cores = num_cores, function(x){ # ri
          mean(exp_table[x, gene])
        }))
        
        mut_perms <- mclapply(1:num_permutations_negative_correlation, mc.cores = num_cores, function(x){ # returns table with rows = vertices, columns = permutations
          perm <- mut_bin[sample(nrow(mut_bin)), ]
          sapply(graph, function(y) mean(perm[y, gene]))
        })
        
        mut_perms <- do.call('cbind', mut_perms)
        
        mut_freq_norm <- mut_freq / sum(mut_freq)
        mut_perms_norm <- sweep(mut_perms, 2, colSums(mut_perms), `/`)
        mut_perms_norm[which(is.na(mut_perms_norm))] <- 0
        avg_exp_norm <- avg_exp / sum(avg_exp)
        
        jsdiv <- suppressMessages(JSD(rbind(avg_exp_norm, mut_freq_norm), est.prob = 'empirical'))
        jsdiv_perm <- suppressMessages(JSD(rbind(avg_exp_norm, t(mut_perms_norm))))
        
        p <- sum(tail(jsdiv_perm[ , 1], -1) < jsdiv) / num_permutations
        
        neg_cor_list[[count]] <- rbind(neg_cor_list[[count]], data.frame(gene, p))
      }
      
      q <- p.adjust(neg_cor_list[[count]]$p, method = 'BH')
      neg_cor_list[[count]] <- cbind(neg_cor_list[[count]], q)
      
    }
  }
  
  ######## COMPUTING GENE SCORES ########
  
  nonsyn_features <- t(nonsyn_mat[,LGG_drivers]) %>% as.data.frame ##CHANGE
  num_complexes <- length(nerve_complexes)
  gene_scores <- vector('list', num_complexes)
  
  for(i in seq(1, num_complexes, 1)){
    message(paste("Assessing genes in complex", i, "of", num_complexes))
    gene_scores[[i]] <- rayleigh_selection(nerve_complexes[[i]], nonsyn_features, num_perm = num_permutations_gene_localization, num_cores = num_cores, one_forms = FALSE, seed = seed)
    if (i == 18){
      sapply(gene_scores[i], '[', 1)
    }
    gene_scores[[i]] <- data.frame("Gene"= rownames(nonsyn_features), "R" = sapply(gene_scores[i], '[', 1), "p" = sapply(gene_scores[i], '[', 2), "q" = sapply(gene_scores[i], '[', 3))
  }
  
  TDAmut_object@gene_scores <- gene_scores
  TDAmut_object@negative_correlations <- neg_cor_list
  
  return(TDAmut_object)
  
}