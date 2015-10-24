require(sna)
require(grid)
require(ggplot2)
require(ggnet)

#rgraph(6, mode = "graph", tprob = 0.5)
net = read.table("net.csv",sep=",",header=TRUE)
net<-net[,-1]
colnames(net)<-rownames(net)<-toupper(letters[1:ncol(net)])
net  




#net = network(net, directed = FALSE)
label<-sample(seq(from = 0, to = 1, by = 0.1), size = 6)
node.size<-rep(20,6)
node.color<-c("yellow","red","cyan","green","orange","pink")
#coord1<-gplot.layout.kamadakawai(net,layout.par=NULL)
ggnet2(net, node.size = node.size,node.color=node.color,mode=coord1,label=label,labels.alpha=1) + geom_text(aes(label=row.names(net)),vjust=-1,alpha=1,fontface="bold",size=10) + 
  geom_text(aes(label=label),fontface="bold",size=5)
net

guides(label=TRUE, size = FALSE)

ggnet2(net, size = 10, color = "red",layout.par = list(seed.coord=coord1))
ggnet2(net, mode = x)

layout.par<-list(exp=20)
gplot(net,usearrows = F,displaylabels = T,label.pos=10,vertex.enclose = T)
gplot.layout.kamadakawai(net, layout.par)


net = network(net, directed = FALSE)
ggnet(net)


ggnet2(net, color = "phono", palette = "Set1", edge.color = c("color", "grey50"))
x = gplot.layout.fruchtermanreingold(net, NULL)
x

ggnet(net, color = "phono")
devtools::install_github("briatte/ggnet")
library(ggnet)
