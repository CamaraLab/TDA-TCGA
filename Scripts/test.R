
cl <- makeCluster(detectCores())
clusterExport(cl=cl, varlist=c("m","cnv","cnv1"))



m<-Matrix(0,100,100,sparse=TRUE)
m<-parSapply(cl,seq_along(1:nrow(m)), function(i) {
  sapply(seq_along(1:nrow(m)), function (j){
  
    if (identical(cnv1[,i],cnv1[,j])) {m[j,i]<-1} else m[j,i]<-0
    #if (identical(1,1)) m[i,j]<<-1
  })
  #x<-which(ans==1)
  
})

parSapply(cl,seq_along(1:ncol(cnv)), function(i) {
  sapply(seq_along(1:nrow(cnv)), function (j){
    
     cnv[i,j]<-as.numeric(cnv[i,j])
    #if (identical(1,1)) m[i,j]<<-1
  })
  
})


stopCluster(cl)
vec<-x

vec<-1:576
g<-TRUE
sapply(vec,function(i) if (x[i]==y[i]) FALSE)
        
    
}
)

)
