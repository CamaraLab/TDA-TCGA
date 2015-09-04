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

p_value_dist<-sapply(scan$name,function (x) {
  print(x)
  results<-list.files(pattern=paste0(x,".*_results_final"))
  if (length(results) !=0) {
    results<-read.csv(results,as.is=T)$q_value
    results<-sum(results<=0.1,na.rm=T)
  } else {results<-NA}
})



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

p_value_dist<-p_value_dist[complete.cases(p_value_dist)]

rownames(scan)<-scan$name
scan$names
head(scan)

y<-cbind(scan[names(q_value_dist),],q_value_dist,q_integrated_dist,p_value_dist)
head(y)



ggplot(y, aes(x=Resolution, y=Gain, color=q_value_dist, label=q_value_dist)) + 
  scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
  geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("q value distribution")

ggplot(y, aes(x=Resolution, y=Gain, color=q_value_dist, label=q_integrated_dist)) + 
  scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
  geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("q integrated distribution")


ggplot(y, aes(x=Resolution, y=Gain, color=q_value_dist, label=p_value_dist)) + 
  scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
  geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("p value distribution")

