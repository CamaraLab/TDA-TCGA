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
  
  plot_exp_ggplot <- function(exp_table,nonsyn_mut_binary=NA,emb,features) { #for internal use
    
    
    for(feature in features){
      if(all(is.null(nonsyn_mut_binary))){
        gg_color <- exp_table[,feature]
        #mid <- mean(exp_table[,feature]) #if adjusting to average expression across all samples
      }
      else{
        #LGG_drivers <- c('CIC.23152','PTEN.5728','ZNF292.23036','SYNE1.23345','FUBP1.8880','NF1.4763','EGFR.1956','NOTCH1.4851','ATRX.546','NIPBL.25836','IDH1.3417')
        nonsyn_mut_binary<-ifelse(nonsyn_mut_binary>0,1,0)
        nonsyn_mut_binary <- nonsyn_mut_binary[order(rownames(nonsyn_mut_binary)),]
        gg_color <- nonsyn_mut_binary[,feature]
      }
      
      plot <- 
        ggplot(emb) +
        geom_point(aes(x=V1, y=V2, color=gg_color), size=0.9) +
        scale_color_gradient(low='blue',high='red') +
        #scale_color_gradient2(midpoint=mean(t(mutload_mat)), low='blue', mid='grey', high='red',space='Lab') +
        ggtitle(feature)

      
      print(plot + theme(panel.background = element_blank()))
    }
    
  }
  
  build.mapper = function(dist, umap_emb, num_intervals, percent_overlap) {
    sink()
    sink('NULL')
    mapper2D(dist, umap_emb, num_intervals, percent_overlap)
  }

  # INPUT AND CLEAN ---------------------------------------------------------

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

  exp_table_top <- exp_table[,order(-apply(exp_table,2,var))][,1:var_threshold]
  dist_matrix <- as.matrix(cor.dist(as.matrix(exp_table_top)))

    # MAPPER PLOTS ---------------------------------------------------------

  # Initialize parameters to build several Mapper complexes
  interval_range <- seq(min_interval, max_interval, by = interval_step)
  #percent_range <- seq(min_percent_overlap, max_percent_overlap, by = percent_step)
  percent_range <- c(33,60,71.4,77.8,81.8,84.6,86.7,88.2)

  # Checks to user input...fix
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
    emb <- umap(dist_matrix)$layout %>% as.data.frame
  }
    else if (filter_method == "PCA") {
      emb <- autoplot(prcomp(dist_matrix))$data[,1:2]
      #autoplot(prcomp(dist_matrix))
    }
      else if (filter_method == "KNN") {
        knn_graph <- nng(dist_matrix, k=k)
        knn_dist <- shortest.paths(knn_graph)
        emb <- cmdscale(knn_dist,2) %>% as.data.frame
      }
      else {
        warning("Not a valid filter method")
      }
  
  # Building nerve complexes
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
  
  sink()
  
  return(mut_objects)
}


