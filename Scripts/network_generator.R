require(sna)
require(grid)
require(ggplot2)
require(ggnet)

#net<-rgraph(6, mode = "graph", tprob = 0.5)
#colnames(net)<-rownames(net)<-toupper(letters[1:ncol(net)])
net = read.table("net1.csv",sep=",",header=TRUE,row.names=1)
coord1<-as.matrix(read.csv("coord1.csv",row.names=1))
#write.csv(net,"net.csv")
#coord1 = gplot.layout.fruchtermanreingold(net, NULL)
#write.csv(coord1,"coord1.csv")

mprob<-sample(seq(from = 0, to = 1, by = 0.1), size = 6)
mprob<-c(100,2000,5,300,100,20)
node.size<-c(20,20,20,20,20,20)
node.color<-c("yellow","pink","cyan","green","orange","grey")
#coord1<-gplot.layout.kamadakawai(net,layout.par=NULL)
ggnet2(net, size = node.size,node.color=as.numeric(mprob*2),mode=coord1,label=mprob,label.size=8,legend.size = 12, legend.position = "right") + 
  geom_text(aes(label=row.names(net)),vjust=-1.2,alpha=1,fontface="bold",size=10)
  #geom_text(aes(label=letters[1:6]),fontface="bold",size=5)


p_matrix<-matrix(0,6,6)
for (i in 1:length(mprob)) {
  for (j in 1:length(mprob)) {
    p_matrix[i,j]<-net[i,j]*mprob[i]*mprob[j]
  }
  
}
colnames(p_matrix)<-rownames(p_matrix)<-colnames(net)

write.csv(p_matrix,"p_matrix3.csv")


con=0
for (i in 1:6) {
  for (j in 1:6) {
    con=con+net[i,j]*x[i]*x[j]
  }
}
  
con

