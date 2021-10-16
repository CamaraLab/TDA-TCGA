#' Creates network representations of the expression space of a patient cohort across grid of Mapper parameters
#'
#' @param TDAmut_object TDAmut class with expression data
#' @param var_threshold number of most variable genes to use from input expression table. By default is 4500.
#' @param filter_method filter function to create embedding used by TDAmapper. Includes UMAP, KNN, and PCA. By default is set to UMAP.
#' @param k if KNN filter function is used, specifies number of nearest neighbors. By default is 30
#' @param min_interval minimum number of intervals in codomain of filter function. By default is 10.
#' @param max_interval maximum number of intervals in codomain of filter function
#' @param interval_step step size between intervals in grid of Mapper parameters. By default is 10.
#' @param min_percent_overlap minimum percentage overlap between intervals. By default is 20.
#' @param max_percent_overlap maximum percentage overlap between intervals
#' @param percent_step step size between percentages in grid of Mapper parameters. By default is 10.
#' @param num_bins an integer controlling clustering within the same level set. By default is 10.
#'
#' @return Returns a TDAmut object populated with nerve complexes.
#'
#' @export

compute_complexes <- function(TDAmut_object, var_threshold = 4500, filter_method = 'KNN', k = 30,
                              min_interval = 10, max_interval = 80, interval_step = 10,
                              min_percent_overlap = 20, max_percent_overlap = 90, percent_step = 10,
                              num_bins = 10) {

  ######## APPLYING FILTER FUNCTION ########

  if (is_empty(TDAmut_object@expression_table)){
    stop('Run create_TDAmut_object first to populate object with expression data')
  }

  exp_table <- TDAmut_object@expression_table

  exp_table_top <- exp_table[ , order(-apply(exp_table, 2, var))][ , 1:var_threshold]
  dist_matrix <- as.matrix(bioDist::cor.dist(as.matrix(exp_table_top)))

  if(filter_method == "UMAP") {
    emb <- umap(dist_matrix)$layout %>% as.data.frame
  }
  else if(filter_method == "PCA") {
    emb <- autoplot(prcomp(dist_matrix))$data[ , 1:2]
  }
  else if(filter_method == "KNN") {
    knn_graph <- cccd::nng(dist_matrix, k = k)
    knn_dist <- shortest.paths(knn_graph)
    emb <- cmdscale(knn_dist, 2) %>% as.data.frame
  }
  else {
    stop("Not a valid filter method. Select 'UMAP', 'PCA', or 'KNN'")
  }

  ######## CREATING NERVE COMPLEXES ########

  interval_range <- seq(min_interval, max_interval, by = interval_step)
  percent_range <- seq(min_percent_overlap, max_percent_overlap, by = percent_step)

  count = 0
  num_complexes <- length(interval_range) * length(percent_range)
  nerve_complexes <- vector('list', length = num_complexes)

  for(i in interval_range){
    for(p in percent_range){
      count = count + 1
      message(paste("Creating Nerve Complex", count, "of", num_complexes))

      invisible(capture.output(mapperObj <- mapper2D(dist_matrix, emb, c(i,i), p, num_bins_when_clustering = num_bins)))
      gg <- nerve_complex(mapperObj$points_in_vertex)
      nerve_complexes[[count]] <- gg
    }
  }

  message("Nerve Complexes Created!")

  TDAmut_object@mapper_intervals <- interval_range
  TDAmut_object@mapper_percents <- percent_range
  TDAmut_object@filter_embedding <- emb
  TDAmut_object@nerve_complexes <- nerve_complexes

  return(TDAmut_object)
}

