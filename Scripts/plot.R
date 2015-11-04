setwd("c:/Users/Udi/Downloads/COAD_CUR/Figures/")

library(ggplot2)
library(grid)

results_file<-"scan_summary_genes_coarse.csv"
png_file<-paste0(results_file,".png")
data<-read.csv(results_file)
data$q_0.15[data$gain==1.5]<-0
data$q_0.15[data$gain==2.5 & data$resolution>=70]<-0

# Genes  graph
ggplot(data, aes(x=factor(resolution), y=factor(gain), label=q_0.15,fill=q_0.15)) + 
  scale_fill_gradient2(low = 'maroon',high = 'blue',guide_legend(title = "#Genes",alpha=0.8),breaks=seq(0,12,2)) +
  geom_tile(alpha=0.8) + theme_minimal() + geom_text(size=4,vjust=-0.1) + ggtitle(label="COAD \n significant genes") + 
  xlab("Resolution") + ylab ("Gain") + geom_text(label=paste0("(",y$first_connected_samples,")"),vjust=1.1,size=3) +
  coord_equal() + theme(axis.text=element_text(size=8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold")) + 
  theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) + 
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_equal()  + ggsave(png_file)


#Mutload graph
ggplot(data, aes(x=factor(resolution), y=factor(gain), fill=mutload,label=as.numeric(round(mutload,2)))) + 
  scale_fill_gradient(low = 'red',high = 'blue',guide_legend(title = "p_value",aplha=0.7)) +
  geom_tile(alpha=0.7) + theme_minimal()+geom_text(size=3) + ggtitle(label="COAD \n Mutational load") + 
  xlab("Resolution") + ylab ("Gain") +
  theme(axis.text=element_text(size=8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold"),panel.grid.major = element_blank()) +
  theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) + 
       scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        coord_equal()  + ggsave("mutational_grid_notrescaled_figure.png")

#Number of samples
ggplot(data, aes(x=factor(resolution), y=factor(gain), label=first_connected_samples,fill=first_connected_samples)) + 
  scale_fill_gradient2(low = 'maroon',high = 'blue',guide_legend(title = "#Samples",alpha=0.8)) +
  geom_tile(alpha=0.8) + theme_minimal() + geom_text(size=4) + ggtitle(label="COAD \n First connected samples") + 
  xlab("Resolution") + ylab ("Gain") +
  coord_equal() + theme(axis.text=element_text(size=8),axis.title = element_text(size=12),plot.title = element_text(size=12,face="bold")) + 
  theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) + 
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_equal()  + ggsave("First_connected_samples_grid.png")
