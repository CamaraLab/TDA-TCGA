setwd("C:/Users/Udi/Google Drive/Columbia/LAB/Rabadan/Presentations/Presentation3/Networks/")
require(sna)
require(grid)
require(ggplot2)
require(ggnet)
require(igraph)

#net<-rgraph(6, mode = "graph", tprob = 0.5)
#colnames(net)<-rownames(net)<-toupper(letters[1:ncol(net)])


mut<-read.csv("c:/Users/Udi/SkyDrive/TCGA_CURATED/COAD_CUR/Mutations/COAD_Full_Mutations_binary.csv",row.names=1)
rownames(mut)<-substring(rownames(mut),1,15)


draw_network<-function (adj_mat) {
  print(adj_mat)
  net =as.matrix(read.table(adj_mat,sep=",",header=TRUE,row.names=1))
  node.color<-ifelse(mut[rownames(net),"KRAS.3845"]==0,"grey","red")
  mut_edge<-rownames(mut)[(mut[,"KRAS.3845"]==1)]
  net<-graph.adjacency(net,diag = TRUE)
  edges<-get.edgelist(net)
  vertices_to_remove<-which(!(names(as.list(V(net))) %in% as.character(edges)))
  net<-delete_vertices(net,vertices_to_remove)
  
  if (length(vertices_to_remove)>1)  {
    node.color<-node.color[-vertices_to_remove]
    mut_edge<-mut_edge[-vertices_to_remove]
    
  }
  
  
  edges<-gsub("\\.",replacement = "-",edges)
  
  edge_index<-which(edges[,1] %in% mut_edge & edges[,2] %in% mut_edge)
  edge_color<-rep(rgb(0.9,0.9,0.9),nrow(edges))
  edge_color[edge_index]<-"blue"
 # plot(net,displayisolates=FALSE)
  ggnet2(net, mode="fruchtermanreingold",edge.size=.2,edge.color=edge_color,node.size =1, label.nodes=T, node.color = node.color)+ ggsave(paste0(adj_mat,".svg"))
  
}




sapply (list.files("."), draw_network)
