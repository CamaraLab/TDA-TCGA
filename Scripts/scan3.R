suppressWarnings({
  suppressMessages ({
    library(getopt,quietly = T,warn.conflicts = FALSE)
    library(stringr,quietly = T,warn.conflicts = FALSE)
    library(ggplot2,quietly = T,warn.conflicts = FALSE)
  })
  
})


spec = matrix(c(
  "connectivity","c",1,"character",
  "permutations","p",2,"integer",
  "matrix", "m",1,"character",
  "cores", "q",1,"integer",
  "chunk", "k",2,"integer",
  "mutload", "y",1,"logical",
  "genes","z",1,"logical",
  "test_mode","t",1,"integer"
  
), byrow=TRUE, ncol=4)

#connectivity_script<-"../../../Google Drive/Columbia/LAB/Rabadan/TCGA-TDA/Scripts/connectivity5.R"

output_dir<-"Results"
dir.create(output_dir)


arg<-list(matrix="LUAD.h5",genes=TRUE,connectivity=NULL,permutations=1,cores=4,chunk=25,mutload=FALSE,genes=TRUE,test_mode=3)

arg<-getopt(spec) #Conmment this line for debug mode
if ( is.null(arg$connectivity ) ) {arg$connectivity = "connectivity5.R"}
if ( is.null(arg$permutations ) ) {arg$permutations = 500}
if ( is.null(arg$cores ) ) {arg$cores = 4}
if ( is.null(arg$chunk ) ) {arg$chunk = 25}
if ( is.null(arg$mutload ) ) {arg$mutload = FALSE}
if ( is.null(arg$genes ) ) {arg$genes = FALSE}
if ( is.null(arg$test_mode ) ) {arg$test_mode = 0}

#arg$test_mode=4

#log_file<-"LUAD_Cor_PCA_15_1.5_LUAD.h5-162926-2015-10-12_log.csv"

extract_sample_size<-function(log_file) {
  log_file<-readLines(log_file)
  line_original<-grep(pattern = "Original",log_file)
  line_first_connected<-grep(pattern = "First connected",log_file)
  sample_original<-as.numeric(substring(str_match(log_file[line_original],":\\d*"),2))
  sample_first_connected<-as.numeric(substring(str_match(log_file[line_first_connected],":\\d*"),2))
  return(c(sample_original,sample_first_connected))
}





#Listing all files in the directory (Networks,genes_results and mutload results)
networks<-gsub('.{5}$', '', list.files(pattern=paste0(".json")))
genes_results_files<-list.files(pattern=paste0(".*_genes_results"))
mutload_results_files<-list.files(pattern=paste0(".*_mutload_results"))

# Inferring network files that have not been processed yet
networks_to_process_genes<-sapply(networks,function (x) { sum(grepl(pattern = x,genes_results_files))} )
networks_to_process_genes<-names(networks_to_process_genes[networks_to_process_genes==0])
networks_to_process_mutload<-sapply(networks,function (x) { sum(grepl(pattern = x,mutload_results_files))} )
networks_to_process_mutload<-names(networks_to_process_mutload[networks_to_process_mutload==0])


#scan<-read.csv("dict.csv",as.is=T)
scan<-data.frame(networks=networks,stringsAsFactors = F)

scan$resolution<-as.numeric(sapply(scan$networks, function (x) {
  strsplit(x,"_")[[1]][4]
}))

scan$gain<-as.numeric(sapply(scan$networks, function (x) {
  strsplit(x,"_")[[1]][5]
}))

scan$genes_connectivity<-scan$networks %in% networks_to_process_genes
scan$mutload_connectivity<-scan$networks %in% networks_to_process_mutload

scan$original_samples<-NA
scan$first_connected_samples<-NA

if (arg$test_mode!=0) {
  scan<-scan[1:arg$test_mode,]
}

if (arg$mutload==TRUE) {
	
# Mutational load connectivity
count<-0
for (file in scan$networks[scan$mutload_connectivity]) {
  count<-count+1  
  print("*********************************************")
  print (paste("Mutload Connectivity for Graph:",file,"-",count,"out of",length(scan$mutload_connectivity)))
  run_line<-paste("Rscript", arg$connectivity, "-p",arg$permutations,"-h TRUE -n",file,"-m",arg$matrix,"-q",arg$cores,"-k",arg$chunk)
  system(run_line)
  
}
#Updating scan table with mutload files
#Extracting information from mut_results files
mutload_results_files<-sapply(scan$networks,function (x) {  
  results<-list.files(pattern=paste0("^",x,".*_mutload_results"))
  if (length(results)==0) {results<-0}
  return(results)
})

p_value_mutload<-sapply(mutload_results_files,function (file) {
  if (length(file)!=1 | is.na(file)) {
    results<-NA
  } else results<-round(read.csv(file,as.is=T)$p_value ,3)
})



#Hyper mutations plot

ggplot(scan, aes(x=factor(resolution), y=gain)) + 
  geom_point(size=5,aes(color=p_value_mutload<=0.05)) + geom_text(label=p_value_mutload,vjust=1.6)+theme_bw() + ggtitle("Mutational Load Connectivity") + 
  guides(color = guide_legend(title = paste("mutload <= 0.05"),
                              title.theme = element_text(size=10,angle=0,color="blue"))) +  ggsave(filename = "mutload_grid.png")    


}

if (arg$genes==TRUE) {
	
#Genes connectivity
count<-0  
for (file in scan$networks[scan$genes_connectivity]) {
  count<-count+1  
  print("*********************************************")
  print (paste("Connectivity for Graph:",file,"-",count,"out of",sum(scan$genes_connectivity)))
  #run_line<-paste("Rscript", arg$arg$connectivity, "-p 500 -h TRUE -n",file,"-m",arg$matrix)
  run_line<-paste("Rscript",arg$connectivity,"-p",arg$permutations, "-t 20 -g 100 -n",file,"-m",arg$matrix,"-q",arg$cores,"-k",arg$chunk)
  system(run_line)
  
  print ("Extracting samples count:")
  log_file<-list.files(pattern=paste0("^",file,".*_log"))
  print(log_file)
  if (length(log_file!=0)){
    log_file<-log_file[1] # Taking first log file in case two exists. I t doesnt matter since first connected component is always the same for the same network.
    ans<-extract_sample_size(log_file)
    scan[scan$networks==file,"original_samples"]<-ans[1]
    scan[scan$networks==file,"first_connected_samples"]<-ans[2]
  } 
}

genes_results_files<-
  sapply(scan$networks,function (x) {
    results<-list.files(pattern=paste0("^",x,".*_genes_results"))
    if (length(results)==0) {results<-NA}
    return(results)
  })


extract_value<-function (files,feature,threshold)  {
  #Gets a feature(p_value) and a threshold (0.05) and extract the number of observations across list of files below that threshold for this feature
  ans<-sapply(files,function (file) {
    if (length(file)!=1 | is.na(file)) {
      results<-NA
    } else {
      results<-read.csv(file,as.is=T)[,feature]
      results<-sum(results<=threshold,na.rm=T)
    }
  })
  return (as.numeric(ans))
}



#genes q_value plot
threshold_range<-c(0.1,0.15,0.2)
for (threshold in threshold_range) {
  
  q_value_dist<-extract_value(genes_results_files,"q_value",threshold)
  title<-paste("Genes_results_q_value <=",threshold, "Permutations=",arg$permutations)
  ggplot(scan, aes(x=resolution, y=gain, color=q_value_dist, label=q_value_dist)) + 
    scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
    geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle(title) +
    ggsave(filename = paste0("Genes_results_q_value","_",threshold,".png"))    
}

#p_value_plot
p_value_dist<-extract_value(genes_results_files,"p_value",0.05)
ggplot(scan, aes(x=resolution, y=gain, color=p_value_dist, label=p_value_dist)) + 
  scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
  geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle("Genes_results_p_value<=0.05") +
  ggsave(filename = paste0("Genes_results_p_value_0.05.png"))    


}

#file.remove ("Rplots.pdf")


#Writing scanner summary file:
write.csv(scan,"scanner_summary.csv")

#Number of samples plot:
ggplot(scan, aes(x=resolution, y=gain, color=first_connected_samples, label=first_connected_samples)) + 
  scale_color_gradient2(low = 'white', mid='yellow', high = 'red') +
  geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle(paste("Original number of samples:",scan$original_samples[1])) +
  ggsave(filename = paste0("First_component_samples.png"))  




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

#files_to_tar<-list.files(pattern="SKCM_")
#zip(zipfile = "a.zip",files = files_to_tar)

write.csv(number_of_events(genes_results_files,"q_value",0.2),"number_q_value_0.2.csv")
write.csv(number_of_events(genes_results_files,"q_value",0.15),"number_q_value_0.15.csv")
write.csv(number_of_events(genes_results_files,"q_value",0.1),"number_q_value_0.1.csv")
write.csv(number_of_events(genes_results_files,"p_value",0.05),"number_p_value_0.05.csv")


print ("Removing files to Results Directory")
csv_files<-list.files(pattern = "*.csv")
png_files<-list.files(pattern = "*.png")

files_to_move_to_results<-c(csv_files,png_files)
x<-file.rename(files_to_move_to_results,paste0("Results/",files_to_move_to_results))
#x<-file.rename(files_to_move_to_results,paste0("Results\",files_to_move_to_results))

if (sum(x)==length(files_to_move_to_results)) {
  print ("All results files moved to Results dir, archiving files")
  tar(paste0("Results_",arg$matrix,".tar.gz"),"Results")
} else {
        print ("This files were not moved to Results dir:")
        print (files_to_move_to_results[!x])
}



# Average p_value
#sapply(unique_genes,function (genes) {

## average<-sapply(scan$networks,function (x)
#  {
#    results<-list.files(pattern=paste0("^",x,".*_results_final"))
#    if (length(results) !=0) {
#      results<-read.csv(results,as.is=T)
#      results<-results$p_value[results$Gene_Symbol=="TP53"]
#    } else {results<-NA}
#    results<-as.numeric(results)
#  })


#})

#plot_grid<-function (data1,feature_name,threshold,title) {

#png(filename = paste0("grid_","_",threshold,".png"))
#  filename = paste0(title,"_",feature_name,'_',threshold,".png")
#  data1<-data1
#  ggplot(scan, aes(x=resolution, y=gain, color=data1, label=data1)) + 
##   scale_color_gradient2(low = 'white', mid='cyan', high = 'black') +
#  geom_point(size=5) + theme_bw() + geom_text(vjust=1.6) + ggtitle(title)


#labs (y="Gain",x="Resolution") 
#+ 
# ggsave(filename,width=length(unique(scan$resolution)),height=length(unique(scan$gain)))


#print (filename)

#dev.off()
#}

#connectivity_log<-readLines(paste0(file_prefix,"_log.csv"))
##connectivity_log<-readLines("LUAD_Cor_PCA_15_1.5_LUAD.h5-162926-2015-10-12_log.csv")
##line_original<-grep(pattern = "Original",connectivity_log)
#line_first_connected<-grep(pattern = "First connected",connectivity_log)
#sample_original<-as.numeric(substring(str_match(connectivity_log[line_original],":\\d*"),2))
#sample_first_connected<-as.numeric(substring(str_match(x[line_first_connected],":\\d*"),2))


