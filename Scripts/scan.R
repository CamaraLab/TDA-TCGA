library(ggplot2)
spec = matrix(c(
  "matrix", "m",1,"character"
  
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode

scan<-read.csv("dict.csv",as.is=T)
count<-0
for (file in scan$file) {
  count<-count+1  
  print (paste("Connectivity for Graph:",file,"-",count,"out of",nrow(scan)))
  run_line<-paste("Rscript connectivity5.R -p 25 -t 20 -g 100 -n",file,"-m","SKCM.h5")
  system(run_line)
}

#Extracting information from final_Results files
q_value_dist<-sapply(scan$file,function (x) {
  print(x)
  results<-list.files(pattern=paste0("^",x,".*_results_final"))
  if (length(results) ==1) {
    results<-read.csv(results,as.is=T)$q_value
    results<-sum(results<=0.1,na.rm=T)
  } else {results<-NA}
})

p_value_dist<-sapply(scan$file,function (x) {
  print(x)
  results<-list.files(pattern=paste0("^",x,".*_results_final"))
  if (length(results) !=0) {
    results<-read.csv(results,as.is=T)$q_value
    results<-sum(results<=0.05,na.rm=T)
  } else {results<-NA}
})


#q_value_dist<-q_value_dist[complete.cases(q_value_dist)]
#p_value_dist<-p_value_dist[complete.cases(p_value_dist)]
q_value_dist
p_value_dist

rownames(scan)<-scan$file
scan$file

head(scan)

y<-cbind(scan[names(q_value_dist),],q_value_dist,p_value_dist)
head(y)

write.csv(y,"y.csv")

ggplot(y, aes(x=Resolution, y=Gain, color=q_value_dist, label=q_value_dist)) + 
  scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
  geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("q value<0.1")


ggplot(y, aes(x=Resolution, y=Gain, color=q_value_dist, label=p_value_dist)) + 
  scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
  geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("p value<0.05")



#Calculates sum of a particular variable for a list of genes across scan results
number_of_events<-function(scan,variable,value) {
  genes<-sapply(scan$file,function (x)
  {
    results<-list.files(pattern=paste0("^",x,".*_results_final"))
    if (length(results) !=0) {
      results<-read.csv(results,as.is=T)
      results<-results$Gene_Symbol[results[,variable]<=value]
    } else {results<-NA}
  })

  genes<-unlist(genes)
  genes<-genes[complete.cases(genes)]
  unique_genes<-unique(genes)
  
  genes1<-sapply(unique_genes,function (x) sum(x==genes))
  return(sort(genes1,decreasing = T))
  
}



write.csv(number_of_events(scan,"q_value",0.2),"number_q_value_0.2.csv")
write.csv(number_of_events(scan,"p_value",0.1),"number_p_value_0.1.csv")


# Average p_value
sapply(unique_genes,function (genes) {
  
 average<-sapply(scan$file,function (x)
  {
    results<-list.files(pattern=paste0("^",x,".*_results_final"))
    if (length(results) !=0) {
      results<-read.csv(results,as.is=T)
      results<-results$p_value[results$Gene_Symbol=="TP53"]
    } else {results<-NA}
    results<-as.numeric(results)
  })
  
  
})



sum(genes=="TP53")



