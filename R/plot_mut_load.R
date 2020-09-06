#' Plots histogram of mutational load in log10 scale and heatmap of pvalues indicating mutational load localization
#'
#' @import reshape2
#'
#' @param TDAmut_object object of class TDAmut with expression data, mutation data, nerve complexes, mutational load, and mutational load localization data.
#' @param histogram_breaks number of bins of mutational load histogram. By default is 100.
#' @param pvalue_figs number of values to display after decimal in localization heatmap
#' @param colorbar_low color of colorbar for low-range values. By default is 'red'
#' @param colorbar_mid color of colorbar for mid-range values. By default is 'purple'
#' @param colorbar_high color of colorbar for upper-range values. By default is 'blue'
#'
#' @return Returns mutational load histogram and heatmap of mutational load localization across grid of nerve complexes
#'
#' @export
#'

plot_mut_load <- function(TDAmut_object, histogram_breaks = 100, pvalue_figs = 3,
                          colorbar_low = 'red', colorbar_mid = 'purple', colorbar_high = 'blue') {

  if (is_empty(TDAmut_object@mutational_load) || is_empty(TDAmut_object@mutational_load_localization)){
    stop('Run compute_mut_load first to populate object with mutational load and its localization across nerve complexes')
  }

  ######## PLOTTING HISTOGRAM ########

  mutload <- TDAmut_object@mutational_load
  hist(log10(mutload), breaks = histogram_breaks, plot=TRUE, main="LGG Mutational Load")

  ######## PLOTTING HEATMAP ########

  intervals <- TDAmut_object@mapper_intervals
  percents <- TDAmut_object@mapper_percents
  pvals <- TDAmut_object@mutational_load_localization

  data <- data.frame(Intervals = rep(intervals, each = length(intervals)), Percents = rep(percents, length(percents)), Values = as.numeric(sapply(pvals, '[' , 2)))

  ggplot(data, aes(x = factor(Intervals), y = factor(Percents), fill = Values, label = round(data$Values, pvalue_figs))) +
    geom_tile(alpha = 0.7) + theme_minimal()+ geom_text(size = 3) +
    scale_fill_gradient2(low = colorbar_low, high = colorbar_high, mid = colorbar_mid, midpoint = 0.5, guide_legend(title = "p-value")) +
    ggtitle("Localization across Mapper Complexes") +
    xlab("Number of Intervals") + ylab ("Percent Overlap") +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12), plot.title = element_text(size = 15,face = "bold"), panel.grid.major = element_blank()) +
    theme(panel.border = element_rect(fill=NA, color="grey", size = 0.5)) +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()
}
