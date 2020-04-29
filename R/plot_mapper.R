#' Plots features on nerve complexes and embeddings
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes, mutational load, and mutational load localization data.
#' @param type type of data for given features. Can take values 'mutation', 'mutational_load', or 'expression'. By default is NULL
#' @param features genes of interest in character vector. Only applies if 'mutation' or 'expression' data is selected. By default is NULL
#' @param which_complexes vector of interval and percent combinations specificying which nerve complexes to plot, such as c((int1, percent1), (int2, percent2) ...). By default is 'All'
#' @param colorbar_low color of colorbar for lower-range values. By default is 'blue'
#' @param colorbar_high color of colorbar for upper-range values. By default is 'red'
#' @param include_embedding plot desired features on embedding used to create nerve complexes. By default is FALSE
#' @param seed seed value when plotting nerve complex. By default is 121.
#'
#' @return Returns plot of nerve complex or embedding with desired features
#'
#' @export

plot_mapper = function(TDAmut_object, type = NULL, features = NULL, which_complexes = 'All', colorbar_low = 'blue', colorbar_high = 'red', include_embedding = FALSE, seed = 123) { 

  if (is_empty(TDAmut_object@nerve_complexes)){
    stop('Run compute_complexes first to populate object with nerve complexes')
  }
  if (type == 'mutational_load' && is.null(TDAmut_object@mutational_load)){
    stop('Run compute_mut_load first to populate object with mutational load values')
  }
  
  mapper_intervals <- TDAmut_object@mapper_intervals
  mapper_percents <- TDAmut_object@mapper_percents
  
  ######## SPECIFYING COMPLEXES TO PLOT ########
  
  if (which_complexes == 'All'){
    mapper_objects <- TDAmut_object@nerve_complexes
    param_pairs <- expand.grid(a = mapper_intervals, b = mapper_percents) %>% as.vector #??????????
    message('Plotting all nerve complexes...')
  }
  else if (type(which_complexes) == 'numeric'){
    param_pairs <- which_complexes
    
    for (i in length(which_complexes)){
      
      param_pair <- which_complexes[i]
      interval <- param_pair[1]
      percent <- param_pair[2]
      
      if (!(interval %in% mapper_intervals) || !(percent %in% mapper_percents)) {
        stop('The following interval and percent pair is not valid: ', param_pair)
      }
      
      # Choosing complex from list which corresponds to given parameter combination
      int_ind <- which(interval %in% mapper_intervals)
      per_ind <- which(interval %in% mapper_percents)
      num <- length(mapper_intervals)
      
      which_complexes[i] <- 1 + (int_ind-1)*(num) + (per_ind - 1)
      
    }
    
    mapper_objects <- TDAmut_object@nerve_complexes[which_complexes]
  }
  else{
    message('Specify valid complexes to plot. Options include: 1) 'All' 2) a vector of interval and percent pairs c( c(int1, percent1), c(int2, percent2), ....)')
  }
  
  count = 0
  
  ######## LABELING NODES OF MAPPER COMPLEXES WITH DESIRED FEATURE ########
  
  for (mapper_obj in mapper_objects){
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
        label = 'Mean Expression'
      }
      else if(type == 'mutation'){
        feature_table <- TDAmut_object@nonsyn_mutations
        label = 'Mutational Frequency'
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
            else if(type == 'mutation'){
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
        title(main = (paste(feature, label)))
        title(sub = paste('2D Intervals:', param_pair[count][1], '% Overlap:', param_pair[count][2]))
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
      title(sub = paste('2D Intervals:', param_pair[count][1], '% Overlap:', param_pair[count][2]))
      image.plot(legend.only = T, zlim = range(M), col = gradient)
      
    }
    
    else {
      plot(adj_graph, layout.auto(adj_graph))
    }
  }
 
  ######## LABELING EMBEDDING WITH DESIRED FEATURE ########
  
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
      
      else if(type == 'mutational_load'){
        mutload <- TDAmut_object@mutational_load
        
        gg_color <- mutload
        
        plot <-
          ggplot(emb) +
          geom_point(aes(x = V1, y = V2, color = gg_color), size = 0.9) +
          labs(color = label_emb) + 
          scale_color_gradient(low = 'blue', high = 'red') +
          #scale_color_gradient2(midpoint=mean(t(mutload_mat)), low='blue', mid='grey', high='red',space='Lab') +
          ggtitle('Mutational Load')
        
        print(plot + theme(panel.background = element_blank()))
      }
    }
  }
}
