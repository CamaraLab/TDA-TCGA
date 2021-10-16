#' #' Plots features on nerve complexes and embeddings
#' #'
#' #' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes, mutational load, and mutational load localization data.
#' #' @param type type of data for given features. Can take values 'mutation', 'mutational_load', 'expression', or 'none' to plot complexes with no label. By default is 'none'
#' #' @param features genes of interest in character vector. Only applies if 'mutation' or 'expression' data is selected. By default is NULL
#' #' @param which_complexes vector of interval and percent combinations specificying which nerve complexes to plot, such as c(int1, percent1, int2, percent2 ...). Can plot all complexes by setting to 'All'. By default is NULL which chooses 3 complexes in the middle of the parameter ranges.
#' #' @param colorbar_low color of colorbar for lower-range values. By default is 'blue'
#' #' @param colorbar_high color of colorbar for upper-range values. By default is 'red'
#' #' @param seed seed value when plotting nerve complex. By default is 121.
#' #'
#' #' @return Returns plot of nerve complex or embedding with desired features
#' #'
#' #' @export

plot_mapper = function(TDAmut_object, type = 'none', features = NULL, which_complexes = NULL, colorbar_low = 'blue', colorbar_high = 'red', include_embedding = FALSE, seed = 121) {

  if (is_empty(TDAmut_object@nerve_complexes)){
    stop('Run compute_complexes first to populate object with nerve complexes')
  }

  num_complexes <- length(Glioma_object@nerve_complexes)
  mapper_intervals <- Glioma_object@mapper_intervals
  mapper_percents <- Glioma_object@mapper_percents
  mapper_objects <- NULL

  ######## SPECIFYING COMPLEXES TO PLOT ########

  if (any(which_complexes == 'All')){
    mapper_objects <- TDAmut_object@nerve_complexes
    param_pairs <- expand.grid(a = mapper_intervals, b = mapper_percents) %>% as.vector
    message('Plotting all nerve complexes')
    mapper_objects <- TDAmut_object@nerve_complexes
  }
  else if (is.null(which_complexes)){
    num_complexes <- length(TDAmut_object@nerve_complexes)

    if (num_complexes <= 2){ # should not happen
      mapper_objects <- TDAmut_object@nerve_complexes
    }

    mid <- ceiling(num_complexes / 2)
    mapper_objects <- TDAmut_object@nerve_complexes[(mid - 1):(mid + 1)]

  }
  else if (is.numeric(which_complexes)) {

    count = 1

    for (param_pair in split(which_complexes, ceiling(seq_along(which_complexes)/2))){

      interval <- param_pair[1]
      percent <- param_pair[2]

      if (!(interval %in% mapper_intervals) || !(percent %in% mapper_percents)) {
        stop('The following interval and percent pair is not valid: ', param_pair)
      }

      # Choosing complex from list which corresponds to given parameter combination
      int_ind <- which(mapper_intervals %in% interval)
      per_ind <- which(mapper_percents %in% percent)
      num <- length(mapper_intervals)

      which_complexes[count:(count + 1)] <- 1 + (int_ind - 1)*(num) + (per_ind - 1)
      count = count + 2

      message(paste('Plotting nerve complex with 2D intervals =', interval ,'and percent overlap =', percent))

    }

    mapper_objects <- TDAmut_object@nerve_complexes[unique(which_complexes)]
  }
  else{
    message("Specify valid complexes to plot. Options include: 1) 'All' 2) a vector of interval and percent pairs list( c(int1, percent1), c(int2, percent2), ....)")
  }

  count = 1

  ######## LABELING NODES OF MAPPER COMPLEXES WITH DESIRED FEATURE ########

  for (mapper_obj in mapper_objects){
    if(!(is.null(features))){

      # determining which feature to plot
      if(type == 'expression'){
        feature_table <- TDAmut_object@expression_table
        label = 'Mean Expression'
      }
      else if(type == 'mutation'){
        feature_table <- TDAmut_object@nonsyn_mutations
        label = 'Mutational Frequency'
      }
      else{
        stop('Type of feature specified is invalid')
      }

      # plotting
      for (feature in features){
          plot_skeleton(mapper_obj, feature_table[, feature])
          title(main = (paste(feature, label)))
      }

    }

    else if(type == 'mutational_load') {
      mutload <- TDAmut_object@mutational_load

      plot_skeleton(mapper_obj, mutload)
      title(main = 'Mutational Load')
    }

    else {

      plot_skeleton(mapper_obj)

    }
  }
}
