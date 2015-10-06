spec = matrix(c(
  "file","f",1,"character"
  "matrix", "m",1,"character"
  "hyper","h",2,"logical"
  "regular","g",2,"logical",
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode
if ( is.null(arg$hyper ) ) {arg$hyper= FALSE}
if ( is.null(arg$regular ) ) {arg$regular= FALSE}



library(ggplot2)
scan<-read.csv("dict.csv",as.is=T)

scan$resolution<-sapply(scan$file, function (x) {
  strsplit(x,"_")[[1]][4]
})

scan$gain<-sapply(scan$file, function (x) {
  strsplit(x,"_")[[1]][5]
})

# Mutational load connectivity
if (arg$hyper) {
  count<-0
  for (file in scan$file[1:5]) {
    count<-count+1  
    print("*********************************************")
    print (paste("Mut load Connectivity for Graph:",file,"-",count,"out of",nrow(scan)))
    run_line<-paste("Rscript connectivity5.R -p 500 -h TRUE -n",file,"-m","SKCM.h5")
    system(run_line)
  }
  
  #Extracting information from mut_results files
  
  
  mutload_results_files<-sapply(scan$file,function (x) {
      results<-list.files(pattern=paste0("^",x,".*_mutload_results"))
    })
  
  p_value_mutload<-sapply(mut_results_files,function (file) {
    results<-read.csv(file,as.is=T)$p_value
  })
  
  
  
  #Hyper mutations plot
  
  ggplot(scan, aes(x=resolution, y=gain)) + 
    geom_point(size=5,aes(color=p_value_mutload<=0.05)) + geom_text(label=p_value_mutload,vjust=1.6)+theme_bw() + ggtitle("p value<0.05") + 
    guides(color = guide_legend(title = "Significant mutload",
                                title.theme = element_text(size=10,angle=0,color="blue")))
  
} 



#Genes connectivity

if (arg$regular) {
  count<-0
  for (file in scan$file) {
    count<-count+1  
    print("*********************************************")
    print (paste("Connectivity for Graph:",file,"-",count,"out of",nrow(scan)))
    run_line<-paste("Rscript connectivity5.R -p 500 -t 20 -g 100 -n",file,"-m","SKCM.h5")
    system(run_line)
  }
  
  genes_results_files<-
    sapply(scan$file,function (x) {
      results<-list.files(pattern=paste0("^",x,".*_results_final"))
    })
  
  
  
  #Extracting information from genes_results files
  p_value_dist<-sapply(genes_results_files,function (file) {
    results<-read.csv(file,as.is=T)$p_value
    results<-sum(results<=0.05,na.rm=T)
  })
  
  q_value_dist<-sapply(genes_results_files,function (file) {
    results<-read.csv(file,as.is=T)$q_value
    results<-sum(results<=0.2,na.rm=T)
  })
  
  #q_value plot
  ggplot(scan, aes(x=resolution, y=gain, color=q_value_dist, label=q_value_dist)) + 
    scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
    geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("q value<0.2")
  
  #p_value_plot
  ggplot(scan, aes(x=resolution, y=gain, color=p_value_dist, label=p_value_dist)) + 
    scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
    geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("p value<=0.05")
}



#Calculates sum of a particular feature for a list of genes across scan results
number_of_events<-function(genes_results_files,feature,threshold) {
  
  genes<-sapply(genes_results_files,function (file)
  {
      results<-read.csv(file,as.is=T)
      results<-results$Gene_Symbol[results[,feature]<=threshold]
  })
  
  genes<-unlist(genes)
  genes<-genes[complete.cases(genes)]
  unique_genes<-unique(genes)
  
  genes_events<-sapply(unique_genes,function (x) sum(x==genes))
  return(sort(genes_events,decreasing = T))
  
}



write.csv(number_of_events(genes_results_files,"q_value",0.2),"number_q_value_0.2.csv")
write.csv(number_of_events(genes_results_files,"p_value",0.05),"number_p_value_0.05.csv")


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


y<-scan$file
y
split_Str(y)

d<-1


t<-rbind(scan,sc)

strsplit(y,"_")[[1]][4]
