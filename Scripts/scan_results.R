
setwd("../../../Documents/TCGA-DATA/test/")
guid<-332438
all_samples_size<-207

genes_results_files<-list.files(pattern=paste0(guid,".*_genes_results.csv"))

networks<-gsub('.{48}$', '', genes_results_files)
resolution<-sapply(strsplit(networks,split = "_"),"[[",4)
gain<-sapply(strsplit(networks,split = "_"),"[[",5)

scan<-data.frame(networks=networks,resolution=NA,gain=NA,original_samples=NA,first_connected_samples=NA,edges_num=NA,samples_threshold=NA,above_samples_threshold=NA,above_gscore_threshold=NA,p_0.05=NA,q_0.1=NA,q_0.15=NA,q_0.2=NA,mutload=NA,parsing_time=NA,connectivity_time=NA,stringsAsFactors = F)
scan$original_samples<-all_samples_size
scan$resolution<-resolution
scan$gain<-gain



#Extract first connected samples  
for (file in scan$networks) {
  log_file<-read.csv(list.files(pattern=paste0(file,".+",guid,".+","_log.csv")),as.is=TRUE)
  first_connected_samples<-log_file[grep("First connected",log_file[,1]),]
  first_connected_samples<-strsplit(first_connected_samples,":")[[1]][2]
  scan[scan$networks==file,]$first_connected_samples<-first_connected_samples
  
  
}


head(scan)
#count<-0
#for (file in scan$networks) {
#  count<-count+1
#  print (paste("Analyzing network:",file,"-",count,"out of",nrow(scan)))
  
#}




################ Number of events per gene###########################

#Calculates sum of a particular feature for a list of genes across scan results
number_of_events<-function(genes_results_files,feature,threshold) {
  
  genes<-sapply(genes_results_files,function (file)
  {
    if (!is.na(file)) {
      results<-read.csv(file,as.is=T)
      results<-results$Gene_Symbol[results[,feature]<=threshold]  
    }
    
  })
  
  
  genes<-unlist(genes)
  genes<-genes[complete.cases(genes)]
  unique_genes<-unique(genes)
  genes_events<-sapply(unique_genes,function (x) sum(x==genes))
  return(sort(genes_events,decreasing = T))
  
}


#Generating number of events summary file
q_value_0.1<-number_of_events(genes_results_files,"q_value",0.1)
q_value_0.15<-number_of_events(genes_results_files,"q_value",0.15)
q_value_0.2<-number_of_events(genes_results_files,"q_value",0.2)
p_value_0.05<-number_of_events(genes_results_files,"p_value",0.05)



n<-max(length(q_value_0.1),length(q_value_0.15),length(q_value_0.2),length(p_value_0.05))
length(q_value_0.1)<-n ; length(q_value_0.15) <-n; length(q_value_0.2) <-n; length(p_value_0.05) <-n

events<-data.frame(
  q_value_0.1,
  q_value_0.15,
  q_value_0.2,
  p_value_0.05
)


colnames(events)<-c("q_value_0.1","q_value_0.15","q_value_0.2","p_value_0.05")
write.csv(events,paste0("number_of_events_",guid,".csv"))

for (file in scan$networks) {
  results_file<-read.csv(list.files(pattern=paste0(file,".+",guid,".+","_genes_results.csv")),as.is=TRUE)
  scan[scan$networks==file,]$p_0.05<-sum(results_file[,"p_value"]<=0.05)
  scan[scan$networks==file,]$q_0.1<-sum(results_file[,"q_value"]<=0.1)
  scan[scan$networks==file,]$q_0.15<-sum(results_file[,"q_value"]<=0.15)
  scan[scan$networks==file,]$q_0.2<-sum(results_file[,"q_value"]<=0.2)
  
}

















write.csv(scan,paste0("scan_summary_",guid,".csv"))





