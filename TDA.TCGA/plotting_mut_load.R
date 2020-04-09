plot_mut_load <- function(mut_objects, mut_load, rescale = NULL){
  
  # Localization across grid of Mapper parameters
  intervals = seq(10, 80, by = 10) #CHANGE
  #percents = seq(15,85,by = 10)
  percents <- c(33,60,71.4,77.8,81.8,84.6,86.7,88.2) #CHANGE
  pval_data <- matrix(sapply(pvals_sub,'[',2),8,8,byrow=TRUE) %>% as.data.frame
  colnames(pval_data) <- percents
  dat2 <- cbind(intervals,pval_data) %>% as.matrix
  mode(dat2) = 'numeric'
  dat2 <- dat2 %>% as.data.frame
  data <- melt(dat2,id.vars='intervals',measure.vars = c(as.character(percents)))
  data <- data.frame(Intervals = rep(intervals,each=8), Percents = rep(percents,8), Values=as.numeric(data$value))
  
  ggplot(data, aes(x=factor(Intervals),y=factor(Percents),fill=Values,label=round(data$Values,3))) +
    geom_tile(alpha=0.7) + theme_minimal()+ geom_text(size=3) +
    scale_fill_gradient2(low='red',high='blue',mid='purple',midpoint=0.5,guide_legend(title="p-value")) +
    ggtitle("Localization across Mapper Complexes") +
    xlab("Number of Intervals") + ylab ("Percent Overlap") +
    theme(axis.text=element_text(size= 8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold"),panel.grid.major = element_blank()) +
    theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()
  
  
  # Plotting Mapper Graphs
  set.seed(121)
  adj_graph <- graph.adjacency(mut_objects$adjacency, mode="undirected")
  #adj_graph <- delete.vertices(adj_graph,degree(adj_graph) < 1)
  points_in_vertex <- mut_objects$points_in_vertex
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