plot_mut_load <- function(mut_objects, mut_load_objects, rescale = NULL){
  
  # Without subsampling
  intervals = seq(10, 80, by = 10) # CHANGE to adapt to mut_objects input
  percents = seq(15,85,by = 10)
  pval_data <- matrix(sapply(mut_load_objects,'[',2),8,8,byrow=TRUE) %>% as.data.frame 
  colnames(pval_data) <- percents
  dat2 <- cbind(intervals,pval_data) %>% as.data.frame
  data <- melt(dat2,id.var='intervals',measure.vars = c(as.character(percents)))
  data <- data.frame(Intervals = rep(intervals,each=8), Percents = rep(percents,8), Values=as.numeric(vals))
  
  ggplot(data, aes(x=factor(Intervals),y=factor(Percents),fill=Values,label=round(data$Values,2))) +
    geom_tile(alpha=0.7) + theme_minimal()+ geom_text(size=3) +
    scale_fill_gradient2(low='red',high='blue',mid='purple',midpoint=0.5,guide_legend(title="p-value")) +
    ggtitle("Localization across Mapper Complexes") +
    xlab("Number of Intervals") + ylab ("Percent Overlap") +
    theme(axis.text=element_text(size= 8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold"),panel.grid.major = element_blank()) +
    theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()
  
  # With subsampling
  
  if(!is.null(rescale)){
    intervals = seq(10, 80, by = 10) # CHANGE to adapt to mut_objects input
    percents = seq(15,85,by = 10)
    pval_data <- matrix(sapply(mut_load_objects,'[',2),8,8,byrow=TRUE) %>% as.data.frame 
    colnames(pval_data) <- percents
    dat2 <- cbind(intervals,pval_data) %>% as.data.frame
    data <- melt(dat2,id.var='intervals',measure.vars = c(as.character(percents)))
    data <- data.frame(Intervals = rep(intervals,each=8), Percents = rep(percents,8), Values=as.numeric(vals))
    
    ggplot(data, aes(x=factor(Intervals),y=factor(Percents),fill=Values,label=round(data$Values,2))) +
      geom_tile(alpha=0.7) + theme_minimal()+ geom_text(size=3) +
      scale_fill_gradient2(low='red',high='blue',mid='purple',midpoint=0.5,guide_legend(title="p-value")) +
      ggtitle("Localization across Mapper Complexes") +
      xlab("Number of Intervals") + ylab ("Percent Overlap") +
      theme(axis.text=element_text(size= 8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold"),panel.grid.major = element_blank()) +
      theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) +
      scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      coord_equal()
    
  }
  
}