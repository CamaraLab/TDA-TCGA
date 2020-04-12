#' Plots features on Mapper graph
#'
#' @param mapper_obj nerve complex on which to plot features
#' @param features genes of interest. By default is NA.
#' @param feature_table table with feature information such as mutation or expression table. By default is NA.
#' @param mutational_load vector of mutational load for each sample. By default is NA.
#'
#' @return Returns labeled plot of nerve complex with given feature
#'
#' @export

plot_mapper = function(mapper_obj, features=NA, feature_table=NA, mutational_load=NA) { #mutload of class named numeric

  set.seed(121)
  adj_graph <- graph.adjacency(mapper_obj$adjacency, mode="undirected")
  #adj_graph <- delete.vertices(adj_graph,degree(adj_graph) < 1)
  points_in_vertex <- mapper_obj$points_in_vertex
  V(adj_graph)$size <- log2(as.numeric(lapply(points_in_vertex, length))+1)
  E(adj_graph)$size <- rep(0, gsize(adj_graph))
  V(adj_graph)$label <- ""
  V(adj_graph)$frame.color <- "black"

  if(!all(is.na(features))){

    M <- matrix(0, nrow=length(points_in_vertex), ncol=length(features), dimnames=list(NULL, features))

    for (feature in features){
      for (i in 1:length(points_in_vertex)) {
        points = points_in_vertex[[i]]
        if (length(points) == 1) {
          M[i, ] <- 0
        }
        else {
          M[i, ] <- length(which(feature_table[points,feature] != 0)) / length(points_in_vertex) # mutational frequency
        }
      }
      gradient <- colorRampPalette(c('blue','red'))
      gradient <- gradient(length(points_in_vertex))
      sample_order <- order(M[,1])
      colors <- vector('double', length = length(points_in_vertex))
      colors[sample_order] <- gradient

      V(adj_graph)$color <- colors
      set.seed(123)
      plot(adj_graph, layout.auto(adj_graph))
      title(main = (paste0(feature,' mutational frequency')))
      image.plot(legend.only = T, zlim=range(M), col=gradient)
    }
  }

  else if(!all(is.na(mutational_load))) {
    points_in_vertex <- mapper_obj$points_in_vertex

    M <- matrix(0, nrow=length(points_in_vertex), ncol=1)

    for (i in 1:length(points_in_vertex)) {
      points <- points_in_vertex[[i]]
      M[i, ] <- sum(mutational_load[points])
    }

    gradient <- colorRampPalette(c('blue','red'),bias=1)
    gradient <- gradient(length(points_in_vertex))
    sample_order <- order(M[,1])
    colors <- vector('double', length = length(points_in_vertex))
    colors[sample_order] <- gradient

    V(adj_graph)$color <- colors
    plot(adj_graph, layout.auto(adj_graph))
    image.plot(legend.only = T, zlim=range(M), col=gradient)
  }

  else {
    plot(adj_graph, layout.auto(adj_graph))
  }
}
