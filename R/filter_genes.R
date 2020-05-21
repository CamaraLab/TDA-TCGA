#' Filter genes considered for compute_localization by mutation frequency, fraction of nonsynonymous mutations, and negative correlations between expression and mutation data
#' 
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes
#' @param freq_threshold threshold mutation frequency. Genes below this value are not considered in analysis. By default is 0.02
#' @param top_nonsyn_fraction number of genes to keep with greatest nonsyn/nonsyn+syn fraction. By default is 350.
#' @param negative_correlations if TRUE, assesses negative correlations between mutations rates and expression of a gene. By default is TRUE.
#' @param num_permutations if negative_correlations = TRUE, number of permutations used to create null distribution when identifying negative correlations between expression and mutation data. By default is 2000
#' @param num_cores number of cores to be used in computation. By default is 1.
#' @param seed integer specifying the seed used to initialize the generator of permutations. By default is 121
#' 
#' @import philentropy
#' @import RayleighSelection
#' 
#' @return Returns TDAmut object populated with filtered genes' Laplacian scores, p values, and q values for localization and negative correlations
#' 
#' @export

filter_genes <- function(TDAmut_object, freq_threshold = 0.01, top_nonsyn_fraction = 350, negative_correlations = TRUE, num_permutations = 2000, num_cores = 1, seed = 121) {
  
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
  
  nonsyn_bin <- ifelse(nonsyn_mat > 0, 1, 0) %>% as.data.frame
  freq <- colMeans(nonsyn_bin) %>% as.data.frame
  nonsyn_mat <- nonsyn_mat[ , freq > freq_threshold]
  
  if (ncol(nonsyn_mat) >= top_nonsyn_fraction){
    nonsyn_frac <- colSums(nonsyn_mat) / colSums(nonsyn_mat + syn_mat[ , colnames(nonsyn_mat)])
    nonsyn_mat <- nonsyn_mat[ , order(nonsyn_frac, decreasing = TRUE)][ , 1:top_nonsyn_fraction]
  }
  
  genes <- colnames(nonsyn_mat)
  ###############################
  
  dist4p <- function(m) { # stack exch adapted for 0s in prob dist
    ncol <- ncol(m)
    
    m2 <- ifelse(is.finite(log2(2*m)), m * log2(2 * m), 0)
    
    xlogx <- matrix(colSums(m2), ncol, ncol)
    xlogx2 <- xlogx + t(xlogx)
    xlogx2[upper.tri(xlogx2, diag=TRUE)] <- 0
    
    xx <- lapply(seq_len(ncol)[-1], function(i, m) {
      j <- seq_len(i - 1L)
      xy <- m[, i] + m[, j, drop=FALSE]
      xy <- ifelse(is.finite(log2(xy)), xy*log2(xy),0)
      colSums(xy)
    }, m)
    
    xylogxy <- matrix(0, ncol, ncol)
    xylogxy[upper.tri(xylogxy, diag=FALSE)] <- unlist(xx)
    
    (0.5 * (xlogx2 - t(xylogxy)))
  }
  
  worker <- function(genes){
    
    exp_table_ofinterest <- exp_table[ , genes, drop = FALSE] %>% as.matrix
    mut_bin_ofinterest <- mut_bin[ , genes, drop = FALSE] %>% as.matrix
    
    for (graph in vertices_in_graphs){
      
      count = count + 1
      # message('Graph ', count)
      
      avg_exp <- sapply(graph, function(x) colMeans(exp_table_ofinterest[x, , drop = FALSE], na.rm = TRUE))
      
      if (length(genes) == 1){
        avg_exp_norm <- list(avg_exp / sum(avg_exp))
      } else {
        avg_exp_norm <- sweep(avg_exp, 1, rowSums(avg_exp), `/`) 
        avg_exp_norm <- avg_exp_norm %>% split(seq(nrow(avg_exp_norm)))
      }
      
      p_list <- rep(NA, length(genes)) %>% as.list
      mut_freq <- matrix(NA, num_permutations, length(graph))
      jsdiv_mat <- matrix(NA, num_permutations, num_permutations)
      
      for (gene in seq_along(genes)) {
        
        message('Assessing ', genes[gene])
        
        mut_freq <- pushCpp(as.numeric(mut_bin_ofinterest[ , gene]), graph, num_permutations)
        mut_freq_norm <- sweep(mut_freq[[1]], 1, rowSums(mut_freq[[1]]), `/`) 
        
        avg_exp_norm[[gene]][which(is.na(avg_exp_norm[[gene]]))] <- 0
        
        jsdiv_mat <- dist4p(cbind(avg_exp_norm[[gene]], t(mut_freq_norm)))
        
        jsdiv <- jsdiv_mat[2,1]
        p_list[[gene]] <- sum(tail(jsdiv_mat[ , 1], -2) <= jsdiv) / num_permutations
      }
      
      neg_cor_list[[count]] <- data.frame('Gene' = genes, 'p' = unlist(p_list))
      
      
    }
    return(neg_cor_list)
  }
  
  
  if (!all(genes %in% colnames(exp_table))){
    missing_genes <- genes[!genes %in% colnames(exp_table)]
    genes <- genes[!genes %in% missing_genes]
    message('Not considering the following genes with no expression data: ', paste("'", missing_genes, "'", collapse = ", ", sep = ""))
  }
  
  # Isolating sample indices within vertices of each nerve complex
  vertices_in_graphs <- lapply(nerve_complexes, function(x) {
    x <- x$points_in_vertex
    return(x)
  })
  
  neg_cor_list <- replicate(length(nerve_complexes), data.frame('Gene' = NULL, 'p' = NULL), simplify = F) %>% set_names(1:length(nerve_complexes))
  count = 0
  RNGkind("L'Ecuyer-CMRG") # For reproducible random number generation with forking parallelization
  set.seed(seed)
  
  # Splitting genes among cores
  work <- list()
  
  if (num_cores > length(genes)) { 
    num_cores <- length(genes)
  }
  
  if (num_cores == 1 || length(genes) == 1) {
    work[[1]] <- genes
  } else {
    wv <- floor(length(genes)/num_cores)
    wr <- length(genes) - wv*num_cores
    if (wr>0) {
      for (m in 1:wr) {
        work[[m]] <- (genes[(1+(m-1)*(wv+1)):(m*(wv+1))])
      }
      for (m in (wr+1):num_cores) {
        work[[m]] <- (genes[(1+wr+(m-1)*wv):(wr+m*wv)])
      }
    } 
    else {
      for (m in 1:num_cores) {
        work[[m]] <- (genes[(1+(m-1)*wv):(m*wv)])
      }
    }
  }
  
  core_results <- mclapply(work, worker, mc.cores = num_cores, mc.silent = FALSE)
  # core_results <- lapply(work, worker)
  
  no_q <- do.call(function(...) mapply(bind_rows, ..., SIMPLIFY = FALSE), args = (core_results))
  
  neg_cor_list <- lapply(no_q, function(x) {
    q <- p.adjust(x$p, method = 'BH')
    x <- cbind(x, q)
  })
  
  
  TDAmut_object@filtered_genes <- genes
  TDAmut_object@negative_correlations <- neg_cor_list
  
  return(TDAmut_object)
  
}

