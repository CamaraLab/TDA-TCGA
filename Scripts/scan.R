return_source<-function (x) {
  return(x[[5]])
}
return_lab<-function (x) {
  return(x[[7]])
}


scan<-read.csv("coad_scan.csv",as.is=T)
scan$name<-paste0(scan$Resolution,"_",scan$Gain)
scan$source<-sapply(strsplit(scan$Link,"/"),return_source)
scan$lab<-sapply(strsplit(scan$Link,"/"),return_lab)

##Downloading
#for (i in 1:nrow(scan)) {
##  source<-scan$source[i]
#  lab<-scan$lab[i]
#  name<-scan$name[i]
#  print (paste("Downloading",scan$name[i],"---",i,"out of",nrow(scan)))
#  run_line<-paste("python download.py",source,lab,name)
#  system (run_line)
#}

#Conenctivity
count<<-0
sapply(scan$name, function (x) {
  count<<-count+1 
  print (paste("Connectivity",count,"---","out of",nrow(scan)))
  run_line<-paste("Rscript connectivity5.R -m COAD.h5 -p 500 -t 20 -g 100 -r 2.75 --maf PROCESSED_hgsc.bcm.edu_COAD.IlluminaGA_DNASeq.1.somatic.v.2.1.5.0.maf -n",x)
  system(run_line)
})

#Extracting information from final_Results files
q_value_dist<-sapply(scan$name,function (x) {
  print(x)
  results<-list.files(pattern=paste0(x,".*_results_final"))
  if (length(results) !=0) {
    results<-read.csv(results,as.is=T)$q_value
    results<-sum(results<=0.2,na.rm=T)
  } else {results<-NA}
})

x<-"26_3"
q_integrated_dist<-sapply(scan$name,function (x) {
  print(x)
  results<-list.files(pattern=paste0(x,".*_results_final"))
  if (length(results) !=0) {
    results<-read.csv(results,as.is=T)$q_integrated
    results<-sum(results<=0.2,na.rm = T)
  } else {results<-NA}
})

q_value_dist<-q_value_dist[complete.cases(q_value_dist)]
q_integrated_dist<-q_integrated_dist[complete.cases(q_integrated_dist)]
q_value_dist
q_integrated_dist

rownames(scan)<-scan$name
scan$names
head(scan)

y<-cbind(scan[names(q_value_dist),],q_value_dist)
head(y)

ggplot(y(y$Gain~y$Resolution)
plot(y$Resolution,y$Gain)

?text

text(y$Gain~y$Resolution,labels=y$q_value_dist,pos="4")
?
scatterplot3d(x=y$Resolution,y=y$Gain, z=y$q_value_dist,pch=16, highlight.3d=TRUE,
              type="l", main="3D Scatterplot")







ggplot(y, aes(x=Resolution, y=Gain, color=q_value_dist, label=q_value_dist)) + 
  scale_color_gradient2(low = 'white', mid='cyan', high = 'black', midpoint = 7) +
  geom_point(size=5) + theme_bw() + geom_text(hjust=2)




