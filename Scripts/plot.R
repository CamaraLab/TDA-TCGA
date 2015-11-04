library(ggplot2)
library(grid)
x<-read.csv("scan_summary_genes_coarse.csv")
y<-read.csv("scan_summary_mutload.csv")

#Number of samples per graph plot:
ggplot(x, aes(x=factor(resolution), y=factor(gain), label=q_0.15,fill=q_0.15)) + 
  scale_fill_gradient2(low = 'maroon',high = 'blue',guide_legend(title = "#Genes",label = TRUE)) +
  geom_tile() + theme_minimal() + geom_text(size=7,vjust=-0.1) + ggtitle(label="COAD \n significant genes") + 
  xlab("Resolution") + ylab ("Gain") + geom_text(label=paste0("(",y$first_connected_samples,")"),vjust=1,size=6) +
  coord_equal() + theme(axis.text=element_text(size=10),axis.title = element_text(size=15),plot.title = element_text(size=20,face="bold")) + 
  theme(panel.border=element_rect(fill=NA,color="grey")) + 
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))


y<-read.csv("scan_summary_mutload.csv")

ggplot(y, aes(x=factor(resolution), y=factor(gain), fill=mutload)) + 
  scale_fill_continuous(low = 'red',high = 'blue',guide_legend(title = "p_value")) +
  geom_tile(alpha=0.7) + theme_minimal()+geom_text(label=round(y$mutload,2),size=7) + ggtitle(label="COAD \n Mutational load") + 
  xlab("Resolution") + ylab ("Gain") +
  theme(axis.text=element_text(size=10),axis.title = element_text(size=15),plot.title = element_text(size=20,face="bold"),panel.grid.major = element_blank()) +
  theme(panel.border=element_rect(fill=NA,color="grey")) + 
       scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
  
ggplot(x, aes(x=factor(resolution), y=factor(gain), fill=mutload)) + 
  scale_fill_continuous(low = 'red',high = 'blue',guide_legend(title = "p_value",label = TRUE)) +
  geom_tile() + theme_minimal() + geom_text(label=round(y$mutload,2),size=7,vjust=-0.1) + ggtitle(label="COAD \n Mutational load") + 
  xlab("Resolution") + ylab ("Gain") + geom_text(label=paste0("(",y$first_connected_samples,")"),vjust=1,size=6) +
  coord_equal() + theme(axis.text=element_text(size=10),axis.title = element_text(size=15),plot.title = element_text(size=20,face="bold"))

qplot(data=x,x=resolution,y=gain,geom="tile",fill=q_0.15,label=q_0.15)
x$gain

