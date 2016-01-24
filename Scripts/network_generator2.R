#setwd("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/Presentations/Presentation3/results_COAD_385276_genes_rescaled_3")
setwd("C:/Users/Udi/SkyDrive/TCGA_CURATED/Rips/results_LUAD_379098_genes")
require(sna)
require(grid)
require(ggplot2)
require(ggnet)
require(igraph)
library(rhdf5,quietly = T,warn.conflicts = FALSE)
#net<-rgraph(6, mode = "graph", tprob = 0.5)
#colnames(net)<-rownames(net)<-toupper(letters[1:ncol(net)])


h5file<-"../LUAD.h5"
print (paste0("Loading ",h5file, " file to memory"))
mat_non_syn_bin<-h5read(h5file,"Mutations_Binary")
all_samples<-h5read(h5file,"Mutations_Samples")
all_genes<-h5read(h5file,"Mutations_Genes")

rownames(mat_non_syn_bin)<-all_samples
colnames(mat_non_syn_bin)<-all_genes



#mut<-read.csv("c:/Users/Udi/SkyDrive/TCGA_CURATED/COAD_CUR/Mutations/COAD_Full_Mutations_binary.csv",row.names=1)


mut<-mat_non_syn_bin
#rownames(mut)<-substring(rownames(mut),1,15)




draw_network<-function (adj_mat) {
  print(adj_mat)
  net =as.matrix(read.table(adj_mat,sep=",",header=TRUE,row.names=1))
  node.color<-ifelse(mut[rownames(net),"mut_STK11|6794"]==0,"grey","red")
  mut_edge<-rownames(mut)[(mut[,"mut_STK11|6794"]==1)]
  net<-graph.adjacency(net,diag = TRUE)
  edges<-get.edgelist(net)
  vertices_to_remove<-which(!(names(as.list(V(net))) %in% as.character(edges)))
  net<-delete_vertices(net,vertices_to_remove)
  
  if (length(vertices_to_remove)>0)  {
    node.color<-node.color[-vertices_to_remove]
    mut_edge<-mut_edge[-vertices_to_remove]
    
  }
  
  
  edges<-gsub("\\.",replacement = "-",edges)
  
  edge_index<-which(edges[,1] %in% mut_edge & edges[,2] %in% mut_edge)
  edge_color<-rep(rgb(0.9,0.9,0.9),nrow(edges))
  #edge_color<-rep(col2rgb("grey",0.1),nrow(edges))
  edge_color[edge_index]<-"blue"
  edge_size<-rep(.2,nrow(edges))
  #edge_size[edge_index]<-0.5
  edge_alpha<-rep(1,nrow(edges))
  edge_alpha[edge_index]<-1
  
 # plot(net,displayisolates=FALSE)
  ggnet2(net, mode="fruchtermanreingold",edge.size=edge_size,edge.color=edge_color,node.size =2, label.nodes=T, node.color = node.color,edge.alpha=edge_alpha) + 
   ggsave(paste0(adj_mat,".svg"))
  
}


#a<-"adj_mat_0.862537721898165_379098.csv"
a<-c("adj_mat_0.00166689978570123_379098.csv","adj_mat_0.283277350505545_379098.csv","adj_mat_0.532225066196089_379098.csv","adj_mat_0.862537721898165_379098.csv","adj_mat_1.0034453588796_379098.csv")

sapply (a, draw_network)
#sapply (list.files(".","adj_mat"), draw_network)





