#' Plots features on nerve complexes and embeddings
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes, mutational load, and mutational load localization data.
#' @param type type of data for given features. Can take values 'mutation', 'mutational_load', or 'expression'. By default is NULL
#' @param features genes of interest in character vector. Only applies if 'mutation' or 'expression' data is selected. By default is NULL
#' @param colorbar_low color of colorbar for lower-range values. By default is 'blue'
#' @param colorbar_high color of colorbar for upper-range values. By default is 'red'
#' @param include_embedding plot desired features on embedding used to create nerve complexes. By default is FALSE
#' @param seed seed value when plotting nerve complex. By default is 121.
#'
#' @return Returns plot of nerve complex or embedding with desired features
#'
#' @export

plot_mapper = function(TDAmut_object, type = NULL, features = NULL, colobar_low = 'blue', colorbar_high = 'red', include_embedding = FALSE, seed = 123) { 

  mapper_obj <- TDAmut_object@nerve_complexes
  
  # CHOOSE WHICH MAPPER OBJ TO MOVE FORWARD WITH #
  
  adj_graph <- graph.adjacency(mapper_obj$adjacency, mode="undirected")
  #adj_graph <- delete.vertices(adj_graph,degree(adj_graph) < 1)
  points_in_vertex <- mapper_obj$points_in_vertex
  V(adj_graph)$size <- log2(as.numeric(lapply(points_in_vertex, length))+1)
  E(adj_graph)$size <- rep(0, gsize(adj_graph))
  V(adj_graph)$label <- ""
  V(adj_graph)$frame.color <- "black"

  if(!(is.null(features))){
    
    if(type == 'expression'){
      feature_table <- TDAmut_object@expression_table
      label = ' Mean log2(1+TPM)'
    }
    else if(type='mutation'){
      feature_table <- TDAmut_object@nonsyn_mutations
      label = ' Mutational Frequency'
    }
    else{
      stop('Type of data specified is invalid')
    }

    M <- matrix(0, nrow=length(points_in_vertex), ncol=length(features), dimnames=list(NULL, features))

    for (feature in features){
      for (i in 1:length(points_in_vertex)) {
        points <- points_in_vertex[[i]]
        if (length(points) == 1) {
          #M[i, ] <- 0
          M[i, ] <- feature_table[points,feature] # for testing
        }
        else {
          if(type == 'expression'){
            M[i, ] <- mean(feature_table[points,feature])
          }
          else if(type='mutation'){
            M[i, ] <- length(which(feature_table[points,feature] != 0)) / length(points_in_vertex) # mutational frequency
          }
        }
      }
      
      gradient <- colorRampPalette(c(colorbar_low, colorbar_high))
      gradient <- gradient(length(points_in_vertex))
      sample_order <- order(M[,1])
      colors <- vector('double', length = length(points_in_vertex))
      colors[sample_order] <- gradient

      V(adj_graph)$color <- colors
      set.seed(seed)
      plot(adj_graph, layout.auto(adj_graph))
      title(main = (paste0(feature, label)))
      image.plot(legend.only = T, zlim = range(M), col = gradient)
    }
  }

  else if(type == 'mutational_load') {
    mutload <- TDAmut_object@mutational_load
    
    M <- matrix(0, nrow = length(points_in_vertex), ncol = 1)

    for (i in 1:length(points_in_vertex)) {
      points <- points_in_vertex[[i]]
      M[i, ] <- sum(mutload[points])
    }

    gradient <- colorRampPalette(c(colorbar_low, colorbar_high))
    gradient <- gradient(length(points_in_vertex))
    sample_order <- order(M[,1])
    colors <- vector('double', length = length(points_in_vertex))
    colors[sample_order] <- gradient

    V(adj_graph)$color <- colors
    set.seed(seed)
    plot(adj_graph, layout.auto(adj_graph))
    title(main = 'Mutational Load')
    image.plot(legend.only = T, zlim = range(M), col = gradient)
    
  }

  else {
    plot(adj_graph, layout.auto(adj_graph))
  }
 
  if(include_embedding == TRUE){
    
    emb <- TDAmut_object@filter_embedding
    
    if(!is.null(features)){
      
      if(type == 'expression'){
        feature_table <- TDAmut_object@expression_table
        label_emb <- 'log2(1+TPM)'
      }
      else if(type == 'mutation'){
        feature_table <- TDAmut_object@nonsyn_mutations
        feature_table <- ifelse(feature_table > 0, 1, 0)
        label_emb <- 'Mutation (binary)'
      }
      
      for(feature in features){
        gg_color <- feature_table[ , feature]
      
        plot <-
          ggplot(emb) +
          geom_point(aes(x = V1, y = V2, color = gg_color), size = 0.9) +
          labs(color = label_emb) + 
          scale_color_gradient(low = 'blue', high = 'red') +
          #scale_color_gradient2(midpoint=mean(t(mutload_mat)), low='blue', mid='grey', high='red',space='Lab') +
          ggtitle(feature)
        
        print(plot + theme(panel.background = element_blank()))
      }
  }
  
}
