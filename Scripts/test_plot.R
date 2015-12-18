data_summary <- function(x) {
  mu <- mean(x)
  sigma1 <- mu-sd(x)
  sigma2 <- mu+sd(x)
  return(c(y=mu,ymin=sigma1,ymax=sigma2))
}


genes<-c("exp_EGFR|1956","exp_KEAP1|9817","exp_KRAS|3845")
t<-as.data.frame(mat_tpm[,genes])

ggplot(melt(t),aes(variable,value)) + geom_violin(aes(fill=variable)) + stat_summary(fun.data=data_summary,)



x<-read.csv("../COAD_CUR/Results/Results_COAD_342285_genes_fine/number_of_events_corrected_COAD.csv")
x<-x[1:10,]
x$percentage<-round(x$q_value_0.15/49*100)
x$status<-c(1,1,1,2,3)
colnames(x)[1]<-"gene"
x$gene<-factor(x$gene,levels=x$gene)


ggplot(x,aes(gene,percentage,group=1)) + geom_step() + 
  theme(panel.background = element_rect(fill = 'cyan', size = 2)) + 
  theme(axis.text.x = element_text(angle =45,vjust = 0.8,hjust=1,colour=x$status)) +
  ggplot2.custo 
  geom_rect(data=NULL,aes(xmin=Inf,xmax=-Inf,ymin=-Inf,ymax=Inf),fill="cyan") +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax="TP53",ymin=-Inf,ymax=Inf),fill="lightgreen",alpha=0.05)
  



class(x)
head(x)
