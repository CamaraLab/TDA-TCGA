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

#Listing all files in the directory (Networks,genes_results and mutload results)
networks<-gsub('.{5}$', '', list.files(pattern=paste0(".json")))
genes_results_files<-list.files(pattern=paste0(".*_genes_results"))
mutload_results_files<-list.files(pattern=paste0(".*_mutload_results"))

networks_to_process_genes<-sapply(networks,function (x) { sum(grepl(pattern = x,genes_results_files))} )
networks_to_process_genes<-names(networks_to_process_genes[networks_to_process_genes==0])
networks_to_process_mutload<-sapply(networks,function (x) { sum(grepl(pattern = x,mutload_results_files))} )
networks_to_process_mutload<-names(networks_to_process_mutload[networks_to_process_mutload==0])



#scan<-read.csv("dict.csv",as.is=T)

scan<-data.frame(networks=networks,stringsAsFactors = F)
#scan$networks<-gsub('.{5}$', '', scan$networks)

scan$resolution<-as.numeric(sapply(scan$networks, function (x) {
  strsplit(x,"_")[[1]][4]
}))

scan$gain<-as.numeric(sapply(scan$networks, function (x) {
  strsplit(x,"_")[[1]][5]
}))

scan$genes_connectivity<-scan$networks %in% networks_to_process_genes
scan$mutload_connectivity<-scan$networks %in% networks_to_process_mutload

# Mutational load connectivity
count<-0
for (file in scan$networks[scan$mutload_connectivity]) {
    count<-count+1  
    print("*********************************************")
    print (paste("Mut load Connectivity for Graph:",file,"-",count,"out of",length(scan$mutload_connectivity)))
    run_line<-paste("Rscript connectivity5.R -p 500 -h TRUE -n",file,"-m","SKCM.h5")
    system(run_line)
  }
#Updating scan table with mutload files
#Extracting information from mut_results files
mutload_results_files<-sapply(scan$networks,function (x) {  
    results<-list.files(pattern=paste0("^",x,".*_mutload_results"))
    if (length(results)==0) {results<-NA}
    return(results)
    })
  
p_value_mutload<-sapply(mutload_results_files,function (file) {
    if (length(file)!=1 | is.na(file)) {
      results<-NA
    } else results<-read.csv(file,as.is=T)$p_value 
  })
  
  
  
#Hyper mutations plot
  png('mutload_grid.jpg')
  ggplot(scan, aes(x=factor(resolution), y=gain)) + 
    geom_point(size=5,aes(color=p_value_mutload<=0.05)) + geom_text(label=p_value_mutload,vjust=1.6)+theme_bw() + ggtitle("Mutational Load Connectivity") + 
    guides(color = guide_legend(title = "mutload<=0.05",
                                title.theme = element_text(size=10,angle=0,color="blue")))

  dev.off()



#Genes connectivity
count<-0  
for (file in scan$networks[scan$genes_connectivity]) {
    count<-count+1  
    print("*********************************************")
    print (paste("Connectivity for Graph:",file,"-",count,"out of",sum(scan$genes_connectivity)))
    run_line<-paste("Rscript connectivity5.R -p 500 -t 20 -g 100 -n",file,"-m","SKCM.h5")
    system(run_line)
  }
  
  genes_results_files<-
    sapply(scan$networks,function (x) {
      results<-list.files(pattern=paste0("^",x,".*_genes_results"))
      if (length(results)==0) {results<-NA}
      return(results)
    })

  #Extracting information from genes_results files
  p_value_dist<-sapply(genes_results_files,function (file) {
    print (length(file))
    if (length(file)!=1 | is.na(file)) {
      results<-NA
    } else {
      results<-read.csv(file,as.is=T)$p_value
      results<-sum(results<=0.05,na.rm=T)
    }
  })
  
  q_value_dist<-sapply(genes_results_files,function (file) {
    if (length(file)!=1 | is.na(file)) {
      results<-NA
    } else {
      results<-read.csv(file,as.is=T)$q_value
      results<-sum(results<=0.2,na.rm=T)
    }
  })
  
  #q_value plot
  ggplot(scan, aes(x=factor(resolution), y=gain, color=q_value_dist, label=q_value_dist)) + 
    scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
    geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("q value<=0.2")
  
  #p_value_plot
  ggplot(scan, aes(x=factor(resolution), y=gain, color=p_value_dist, label=p_value_dist)) + 
    scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
    geom_point(size=5) + theme_bw() + geom_text(hjust=2) + ggtitle("p value<=0.05")




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
  
 average<-sapply(scan$networks,function (x)
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


y<-scan$networks
y
split_Str(y)

d<-1


t<-rbind(scan,sc)

strsplit(y,"_")[[1]][4]
