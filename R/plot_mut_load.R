#' Plots histogram of mutational load in log10 scale and heatmap of pvalues indicating mutational load localization
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
                          heatmap_color_low = 'red', heatmap_color_mid = 'purple', heatmap_color_high = 'blue') {
  
  ######## HISTOGRAM ########
  mutload <- TDAmut_object@mutational_load
  hist(log10(mutload),breaks=100,plot=TRUE,main="LGG Mutational Load")
  
  ######## HEATMAP ########
  intervals <- TDAmut_object@mapper_intervals
  percents <- TDAmut_object@mapper_percents
  pvals <- TDAmut_object@mutational_load_localization
  
  # throwaway
  intervals <- seq(10, 80, by = 10)
  #percents <- seq(15,85,by = 10)
  percents <- c(33,60,71.4,77.8,81.8,84.6,86.7,88.2)
  
  pval_data <- matrix(sapply(pvals,'[',2),length(intervals),length(percents),byrow=TRUE) %>% as.data.frame
  colnames(pval_data) <- percents
  dat2 <- cbind(intervals,pval_data) %>% as.matrix
  mode(dat2) = 'numeric'
  dat2 <- dat2 %>% as.data.frame
  data <- melt(dat2,id.vars = 'intervals',measure.vars = c(as.character(percents)))
  data <- data.frame(Intervals = rep(intervals, each = length(intervals)), Percents = rep(percents, each = length(percents)), Values = as.numeric(data$value))
  
  ggplot(data, aes(x=factor(Intervals),y=factor(Percents),fill=Values,label=round(data$Values, pvalue_figs))) +
    geom_tile(alpha=0.7) + theme_minimal()+ geom_text(size=3) +
    scale_fill_gradient2(low = colorbar_low, high = colorbar_high, mid = colorbar_mid, midpoint = 0.5, guide_legend(title="p-value")) +
    ggtitle("Localization across Mapper Complexes") +
    xlab("Number of Intervals") + ylab ("Percent Overlap") +
    theme(axis.text=element_text(size= 8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold"),panel.grid.major = element_blank()) +
    theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()
}