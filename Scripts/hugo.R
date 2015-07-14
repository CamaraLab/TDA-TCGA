#hugo<-read.delim("../Downloads/hugo_genes.txt",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
#hugo2<-read.delim("../../../Downloads/hugo_genes2.txt",header = TRUE,row.names = 1,stringsAsFactors = FALSE)

#colnames(hugo)<-c("prev","Synonyms","Entrez","EntrezNCBI")
#head(hugo)

x<-data.frame(prev=NULL,Synonyms=NULL,Entrez=NULL,EntrezNCBI=NULL)
y<-NULL

y<-apply(hugo,1, function (x){
  if (identical(hugo[i,3],hugo[i,4])) id<-hugo[3] 
  
})

for (i in 1:nrow(hugo)) 
  if (identical(hugo[i,3],hugo[i,4])) y<-hugo[i,3] 
if (!identical(hugo[i,3],hugo[i,4])) {
  if (!is.na(hugo[i,3]) | !is.na(hugo[i,4])) {y<-"?"} else {if (is.na(hugo[i,3]) y<-hugo[i,4]) y<-}
}  

write.csv(x,"not_identical.csv")
View(x)
dim(x)

Anno<-data.frame(Symbol=NULL,Entrez=NULL)
Anno<-cbind(rownames(hugo),hugo$Entrez)
dim(Anno)
head(Anno)

for (i in 1:nrow(Anno))
{
  
  t1<-strsplit(hugo[i,1],", ")
  t1<-as.character(t1[[1]])
  t2<-strsplit(hugo[i,2],", ")
  t2<-as.character(t2[[1]])
  t3<-c(t1,t2)
  if(length(t3)!=0) {
    t4<-cbind(t3,rep(hugo[i,3],length(t3)))
    Anno<-rbind(Anno,t4)
  }
}
rownames(Anno)<-Anno[,1]
