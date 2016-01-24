require("reshape")
require("ggplot2")

delete_background <- theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_blank(),
                           legend.position="none")

#x<-read.csv("c:/Users/Udi/SkyDrive/TCGA_CURATED/STADTRIMMED/Networks/STADTRIM_Extra_Fine_NETWORKS_45_75_4.5_7.5/results_STADTRIM_347546_genes_rescaled_2.6/jsd_value_matrix_347546.csv",row.names=1)
x<-read.csv("c:/Users/Udi/SkyDrive/TCGA_CURATED/COAD_CUR/Results/JSD_results_COAD_341959_genes_rescaled_3_fine/jsd_q_value_matrix_341959.csv",row.names=1)
#xq<-apply(x,2,function(x) p.adjust(x,method="fdr"))
#xq<-xq[g,]
#x<-x[g,]

setwd("c:/Users/Udi/SkyDrive/TCGA_CURATED/COAD_CUR/Results/JSD_results_COAD_341959_genes_rescaled_3_fine/")
write.csv(x,"jsd_value_matrix_341959-1.csv")



rownames(x)<-sapply(strsplit(rownames(x),"|",fixed = TRUE),"[[",1)
xx<-melt(t(x))
xx$X2<-factor(xx$X2,levels=xx$X2)
px<-ggplot(xx,aes(X2,value)) + geom_boxplot(aes(fill=X2),alpha=0.7) + theme(axis.text.x = element_text(angle =45,vjust = 0.8,hjust=1)) + ggtitle ("js_distance")
px + delete_background


y<-read.csv("c:/Users/Udi/SkyDrive/TCGA_CURATED/STADTRIMMED/Networks/STADTRIM_Extra_Fine_NETWORKS_45_75_4.5_7.5/results_STADTRIM_347546_genes_rescaled_2.6/jsd_p_value_matrix_347546.csv",row.names=1)
rownames(y)<-sapply(strsplit(rownames(y),"|",fixed = TRUE),"[[",1)

yy<-melt(t(y))
yy$X2<-factor(yy$X2,levels=yy$X2)
py<-ggplot(yy,aes(X2,value)) + geom_boxplot(aes(fill=X2),alpha=0.7) + theme(axis.text.x = element_text(angle =45,vjust = 0.8,hjust=1)) + ggtitle("js p_value") + geom_hline(yintercept=0.05)
py + delete_background
# + ggsave("jsd_q_value.svg")


yy$value<-p.adjust(yy$value,method = "fdr")
yy$X2<-factor(yy$X2,levels=yy$X2)
py<-ggplot(yy,aes(X2,value)) + geom_boxplot(aes(fill=X2),alpha=0.7) + theme(axis.text.x = element_text(angle =45,vjust = 0.8,hjust=1)) + ggtitle("js q_value") + geom_hline(yintercept=0.15,color="red")
py + theme(axis.text.y = element_text(family="Arial",size=15)) + delete_background + ylim(0,1)
  ggsave("c:/Users/Udi/SkyDrive/TCGA_CURATED/STADTRIMMED/Networks/STADTRIM_Extra_Fine_NETWORKS_45_75_4.5_7.5/jsd_q_value.svg",width=18,height = 8,units="cm")


,face="bold",family="Times")



count<-0
y<-apply(x,2,function(network_genes) {
  count<-count+1
  print(count)
  NA_loc<-which(is.na(network_genes))
  temp<-network_genes
  temp[NA_loc]<-(2)
  #length(temp)
  ans<-sapply(temp,function(x) sum(temp<=x))
  ans[NA_loc]<-NA
  ans<-ans/(length(ans)-length(NA_loc))
  return(ans)
})

write.csv(y,"c:/Users/Udi/SkyDrive/TCGA_CURATED/COAD_CUR//Networks//COAD_Networks_Fine/results_COAD_374366_genes_rescaled_3/jsd_rank.csv")

g<-c("SOX9|6662","APC|324","PIK3CA|5290","ARHGAP5|394","ARFGEF1|10565","TP53|7157","KRAS|3845","VPS13B|157680","SMAD4|4089","TCF7L2|6934","ESRRA|2101")
g<-c(g,"RNF43|54894","KMT2C|58508","FLT3|2322","NEFH|4744","CCDC141|285025","PIK3R1|5295","MYH3|4621","STK11|6794","NCOR1|9611")
xx<-x[g,]
xxx<-melt(t(xx))
yy<-y[g,]

yyy<-melt(t(yy))


xxx$X2<-as.character(xxx$X2)
xxx$X2<-sapply(strsplit(xxx$X2,"|",fixed = TRUE),"[[",1)
xxx$X2<-factor(xxx$X2,levels=xxx$X2)

yyy$X2<-as.character(yyy$X2)
yyy$X2<-sapply(strsplit(yyy$X2,"|",fixed = TRUE),"[[",1)
yyy$X2<-factor(yyy$X2,levels=yyy$X2)

px<-ggplot(xxx,aes(X2,value,colorgroup=1)) + geom_violin() + theme(axis.text.x = element_text(angle =45,vjust = 0.8,hjust=1)) + ggtitle ("js_distance")
px

py<-ggplot(yyy,aes(X2,value,color=X2)) + geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust = 0.8,hjust=1)) + ggtitle("js rank")
py


z<-read.csv("c:/Users/Udi/SkyDrive/TCGA_CURATED/COAD_CUR/Mutations/PROCESSED_MAF_COAD_2015-10-27.maf")
maf<-read.delim("c:/Users/Udi/SkyDrive/TCGA_CURATED/COAD_CUR/Mutations/PROCESSED_MAF_COAD_2015-10-27.maf",header = TRUE,as.is=T,comment.char = "#",sep="\t")

