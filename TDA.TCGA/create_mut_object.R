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
library(Seurat)
library(cccd)
library(maftools)

create_mut_object <- function(exp_table, mut_table, filter_method = UMAP, k = 30,
                              min_interval = 10, max_interval, interval_step = 10,
                              min_percent_overlap = 20, max_percent_overlap, percent_step = 10,
                              var_threshold = 4500, num_cores=1) {


  # FUNCTIONS ---------------------------------------------------------------

  plot.umap = function(emb) {
    
    ggplot(emb) +
      geom_point(aes(x=V1, y=V2, color=0), size=0.9)
      
  }
  
  plot_exp <- function(emb,features,exp_table,threshold){
    
    for(feature in features){
      exp_table <- exp_table[order(row.names(exp_table)), ]
      emb <- emb[order(row.names(emb)), ]
      max <- threshold * max(exp_table[,feature])
      expressed <- row.names(exp_table)[exp_table[,feature] > max]
      nonexpressed <- row.names(exp_table)[exp_table[,feature] < max]  
      
      xpts <- emb[expressed,1]
      x_non <- emb[nonexpressed,1]
      ypts <- emb[expressed,2]
      y_non <- emb[nonexpressed,2]
      
      par(xpd = TRUE, mar=c(1,2,1,10))
      xylim = range(emb)
      xylim = xylim + ((xylim[2]-xylim[1])*0.1)*c(-0.5, 0.5)
      plot(xylim, xylim, type="n", axes=F, frame=F,xlab="",ylab="",asp=1)
      rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
      points(xpts,ypts,pch=16,col='red')
      points(x_non,y_non,pch=16,col='grey')
      title(main=feature)
      legend('topright', inset=c(-.7,0),
             legend=c(paste0("> ",threshold*100,"% Expression"),"< Threshold"),
             col=c('red','grey'),pch=c(16,16), bty='n',x.intersp=0.5)
    }
  }
  
  plot_exp_ggplot <- function(exp_table,nonsyn_mut_binary=NA,emb,features) { #for internal use
    
    
    for(feature in features){
      if(is.na(nonsyn_mut_binary)){
        gg_color <- exp_table[,feature]
        #mid <- mean(exp_table[,feature]) #if adjusting to average expression across all samples
      }
      else{
        #gg_color <- (nonsyn_mut_binary[,feature])
        #test[which(nonsyn_mut_binary[,feature] != 0),feature] <- 1
        nonsyn_mut_binary<-ifelse(nonsyn_mut_binary>0,1,0)
        gg_color <- test[,feature]
        
      }
      
      plot <- 
        ggplot(emb) +
        geom_point(aes(x=V1, y=V2, color=gg_color), size=0.9) +
        scale_color_gradient(low='blue',high='red') +
        #scale_color_gradient2(midpoint=mid, low='blue', mid='grey', high='red',space='Lab') +
        ggtitle(feature)
      
      
      print(plot + theme(panel.background = element_blank()))
    }
    
  }
  
  build.mapper = function(dist, umap_emb, num_intervals, percent_overlap) {
    sink()
    sink('NULL')
    sink()
    mapper2D(dist, umap_emb, num_intervals, percent_overlap)
  }

  # INPUT AND CLEAN ---------------------------------------------------------

  #exp_table <- (read.csv("C:/Users/Adam/Desktop/Camara Lab/TDA_TCGA_Project/Test Data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?"))) #Laptop
  exp_table <- (read.csv("/home/rstudio/documents/Messy_test_data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F)) #Lab
  #exp_table <- (read.csv(exp_table, row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))
  #mut_table <- (read.csv("blah", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))

  if(duplicated(rownames(exp_table))) {
    warning('Cleaning duplicated samples detected in expression data')
    exp_table <- exp_table[!duplicated(rownames(exp_table)),] 
  }
  
  if(duplicated(colnames(exp_table))) {
    warning('Cleaning duplicated genes detected in expression data')
    exp_table <- exp_table[,!duplicated(colnames(exp_table))] 
  }
    
  
  
  # if(all(colnames(exp_table) %in% colnames(mut_table))) { #Gene names can
  #   if(all(colnames(mut_table) %in% colnames(exp_table))) {
  #     message("Matching gene names between expression and mutation tables")
  #   }
  # } else {
  #   print("Unmatched genes between expression and mutation tables")
  # }

  exp_table_top <- exp_table[,order(-apply(exp_table,2,var))][,1:var_threshold]
  dist_matrix <- as.matrix(cor.dist(as.matrix(exp_table_top)))

    # MAPPER PLOTS ---------------------------------------------------------

  # Initialize parameters to build several Mapper complexes
  interval_range <- seq(min_interval, max_interval, by = interval_step)
  percent_range <- seq(min_percent_overlap, max_percent_overlap, by = percent_step)

  # if(!((max_interval - min_interval) %% interval_step)){
  #   interval_range <- c(interval_range, max_interval)
  # }
  # 
  # if(!((max_percent_overlap - min_percent_overlap) %% percent_step)){
  #   percent_range <- c(percent_range, max_percent_overlap)
  # }

  count = 0
  num_complexes <- length(interval_range)*length(percent_range)

  nerveComplexes <- vector('list',length = num_complexes)


  # Filter Functions

  if (filter_method == "UMAP") {

    emb = umap(dist_matrix)$layout %>% as.data.frame

  }

    else if (filter_method == "PCA") {
      emb <- autoplot(prcomp(dist_matrix))$data[,1:2]
      #autoplot(prcomp(dist_matrix))
    }

      else if (filter_method == "KNN") {
        knn_graph <- nng(dist_matrix, k)
        knn_dist <- shortest.paths(knn_graph)
        emb <- cmdscale(knn_dist,2)
      }

      else {
        warning("Not a valid filter method")
      }

  for (i in interval_range){
    for (p in percent_range){
      count = count + 1
      message(paste0("Creating Mapper Complex ", count, " of ", num_complexes))

      mapperObj <- build.mapper(dist_matrix, emb, c(i,i), p)
      gg <- nerve_complex(mapperObj$points_in_vertex)
      nerveComplexes[[count]] <- gg

    }
  }

  mut_objects <- nerveComplexes
  class(mut_objects) <- "DriverMut_TDA"

  return(mut_objects)

  #   # Laplacian eigenmap
  # leim <- LaplacianEigenmaps()
  # lap_emb <- leim@fun(as((exp_table_top), "dimRedData"), leim@stdpars)
  # m_lap <- build.mapper(dist_matrix,list(lap_emb@data@data[,1], lap_emb@data@data[,2]), c(20,20), 50)
  # g_lap <- nerve_complex(m_lap$points_in_vertex)
  # plot_skeleton(g_lap)

}


