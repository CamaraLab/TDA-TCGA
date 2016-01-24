require("ggplot2")
require("getopt")

spec = matrix(c(
  "project_name", "p", 1, "character",
  "guid", "g",1,"integer",
  "sample_size", "s",1,"integer"  
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode

#setwd("C:/Users/Udi/SkyDrive/TCGA_CURATED/LGG_CUR/Results/Results_LGG_312316_genes_coarse/Results_LGG_312316_genes - Copy")
PROJECT_NAME<-arg$project_name

guid<-arg$guid
all_samples_size<-arg$sample_Size

#PROJECT_NAME<-"LGG"
#guid<-2015
#all_samples_size<-512

genes_results_files<-list.files(pattern=paste0(guid,".*_genes_results.csv"))

networks<-gsub('.{40}$', '', genes_results_files)
networks
resolution<-sapply(strsplit(networks,split = "_"),"[[",4)
gain<-sapply(strsplit(networks,split = "_"),"[[",5)

scan<-data.frame(networks=networks,resolution=NA,gain=NA,original_samples=NA,first_connected_samples=NA,edges_num=NA,samples_threshold=NA,above_samples_threshold=NA,above_gscore_threshold=NA,p_0.05=NA,q_0.1=NA,q_0.15=NA,q_0.2=NA,mutload=NA,parsing_time=NA,connectivity_time=NA,stringsAsFactors = F)
scan$original_samples<-all_samples_size
scan$resolution<-resolution
scan$gain<-gain



#Extract first connected samples  
count<-0
for (file in scan$networks) {
  #log_file<-read.csv(list.files(pattern=paste0(file,".+",guid,".+","_log.csv")),as.is=TRUE)
  count<-count+1
  print(count)
  
  print(file)
  log_file<-read.csv(list.files(pattern=paste0(file,".+-",guid,".+","_log.csv")),as.is=TRUE)
  #.+-2015-.+_log.csv
  first_connected_samples<-log_file[grep("First connected",log_file[,1]),]
  first_connected_samples<-strsplit(first_connected_samples,":")[[1]][2]
  scan[scan$networks==file,]$first_connected_samples<-first_connected_samples
    
  
}


for (file in scan$networks) {
  results_file<-read.csv(list.files(pattern=paste0(file,".+",guid,".+","_genes_results.csv")),as.is=TRUE)
  scan[scan$networks==file,]$p_0.05<-sum(results_file[,"p_value"]<=0.05)
  scan[scan$networks==file,]$q_0.1<-sum(results_file[,"q_value"]<=0.1)
  scan[scan$networks==file,]$q_0.15<-sum(results_file[,"q_value"]<=0.15)
  scan[scan$networks==file,]$q_0.2<-sum(results_file[,"q_value"]<=0.2)
  
}






################ Number of events per gene###########################




number_of_events1<-function(genes_results_files1,feature,threshold) {
  results_file_list<-lapply(genes_results_files,function (x) read.csv(x))
  names(results_file_list)<-genes_results_files
  
  genes_freq<-list()
  count<-0
  for (file_name in genes_results_files1) {
    count<-count+1
    #print(count)
    file<-results_file_list[[file_name]]
    new_genes<-setdiff(file$Gene_Symbol,names(genes_freq))
    new_list<-lapply(new_genes,function (x) 0)
    names(new_list)<-new_genes
    genes_freq<-c(genes_freq,new_list)
    for (gene in file$Gene_Symbol) {
      if (file[file$Gene_Symbol==gene,feature]<=threshold) {
        genes_freq[gene]<-genes_freq[[gene]] + 1
      }
    }
  }
  return(genes_freq)
}

  
  

#Generating number of events summary file
q_value_0.1<-number_of_events1(genes_results_files,"q_value",0.1)
q_value_0.15<-number_of_events1(genes_results_files,"q_value",0.15)
q_value_0.2<-number_of_events1(genes_results_files,"q_value",0.2)
p_value_0.05<-number_of_events1(genes_results_files,"p_value",0.05)



events<-t(rbind(
  as.data.frame(q_value_0.1),
  as.data.frame(q_value_0.15),
  as.data.frame(q_value_0.2),
  as.data.frame(p_value_0.05)
))


colnames(events)<-c("q_value_0.1","q_value_0.15","q_value_0.2","p_value_0.05")
write.csv(events,paste0("number_of_events_",guid,".csv"))





write.csv(scan,paste0("scan_summary_",guid,".csv"))

scan$p_0.05[scan$resolution>=40 & scan$gain==1.5]<-0
scan$q_0.1[scan$resolution>=40 & scan$gain==1.5]<-0
scan$q_0.2[scan$resolution>=40 & scan$gain==1.5]<-0
scan$q_0.15[scan$resolution>=40 & scan$gain==1.5]<-0

guid<-352339
setwd("C:/Users/Udi/SkyDrive/TCGA_CURATED/LGG_CUR/Results/Results_LGG_352339_genes_fine/")
scan<-read.csv("scan_summary.csv")

#FIRST CONNECTED  SAMPLES GRID
svg_plot_samples<-paste0("First_connected_samples_grid_",guid,".svg")
plot_samples<-ggplot(scan, aes(x=factor(resolution), y=factor(gain), label=as.numeric(first_connected_samples),fill=as.numeric(first_connected_samples)))
plot_samples<-plot_samples + scale_fill_gradient2(low = 'maroon',high = 'blue',guide_legend(title = "#Samples",alpha=0.8)) +
  geom_tile(alpha=0.8) + theme_minimal() + geom_text(size=4) + ggtitle(label=paste(PROJECT_NAME,"\n First connected samples")) + 
  xlab("Resolution") + ylab ("Gain") +
  coord_equal() + theme(axis.text=element_text(size=8),axis.title = element_text(size=12),plot.title = element_text(size=12,face="bold")) + 
  theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) + 
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_equal()  + ggsave(svg_plot_samples)


q_threshold_range<-c(0.1,0.15,0.2)
for (threshold in q_threshold_range) {
  q_value_dist<-scan[,paste0("q_",threshold)]
  genes_svg_file<-paste0("Genes_results_q_value","_",threshold,"_",guid,".svg")
  ggplot(scan, aes(x=factor(resolution), y=factor(gain), label=q_value_dist,fill=q_value_dist)) + 
    scale_fill_gradient2(low = 'maroon',high = 'blue',guide_legend(title = "#Genes",alpha=0.8),breaks=seq(0,12,2)) +
    geom_tile(alpha=0.8) + theme_minimal() + geom_text(size=4,vjust=-0.1) + ggtitle(label=paste0(PROJECT_NAME," \n significant genes_",threshold,"\n (# of  samples)")) + 
    xlab("Resolution") + ylab ("Gain") + geom_text(label=paste0("(",scan$first_connected_samples,")"),vjust=1.1,size=3) +
    coord_equal() + theme(axis.text=element_text(size=8),axis.title = element_text(size=12),plot.title = element_text(size=15,face="bold")) + 
    theme(panel.border=element_rect(fill=NA,color="grey",size = 0.5)) + 
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()  + ggsave(genes_svg_file)
 
  
}

#connectivity_p_value_plot
p_value_dist<-scan$p_0.05
ggplot(scan, aes(x=resolution, y=gain, label=p_value_dist)) + 
  #scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
  geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle("Genes_results_p_value<=0.05") +
  ggsave(filename = paste0("Genes_results_p_value_0.05_",guid,".svg"))    




