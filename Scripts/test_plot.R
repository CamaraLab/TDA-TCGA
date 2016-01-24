delete_background <- theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_blank(),
                           legend.position="none")

j<-read.csv("../../..//LGG_CUR//Results/Results_LGG_352339_genes_fine/number_of_events_352339.csv",as.is=T)

colnames(j)[1]<-"gene"
jsub<-arrange(j,desc(q_value_0.15))
jsub<-filter(jsub,q_value_0.15>0)
genes<-jsub$gene

networks_num<-length(scan$networks)

jsub$q15_precentage<-round(jsub$q_value_0.15/networks_num*100)
jsub

jsub$gene<-factor(jsub$gene,levels=jsub$gene)

jsub$gene<-as.character(jsub$gene)
p<-ggplot(jsub,aes(gene,q15_precentage,group=1)) + geom_step(col=rgb(1,0.65,0),size=1) + delete_background + 
  theme(panel.background = element_rect(fill = rgb(0.96,0.89,0.84), size = 2)) + 
  theme(axis.text.x = element_text(angle =45,vjust = 0.8,hjust=1,colour=x$status)) + ylim (0,100) + 
  scale_x_discrete(limits=genes)
p




rgb(244,227,215)

class(x)
head(x)
