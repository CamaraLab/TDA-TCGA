require(sna)
require(grid)
require(ggplot2)
require(ggnet)
#net = rgraph(6, mode = "graph", tprob = 0.5)
#net = network(net, directed = FALSE)
label<-colnames(net)<-sample(seq(from = 0, to = 1, by = 0.1), size = 6)
node.size<-rep(c(20,40),3)
#coord1<-gplot.layout.kamadakawai(net,layout.par=NULL)
ggnet2(net, node.size = node.size,mode=coord1,label=label,label.color="red")  
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
