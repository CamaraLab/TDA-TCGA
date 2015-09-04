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

for (i in 1:nrow(scan)) {
  source<-scan$source[i]
  lab<-scan$lab[i]
  name<-scan$name[i]
  print (paste("Downloading",scan$name[i],"---",i,"out of",nrow(scan)))
  run_line<-paste("python download.py",source,lab,name)
  system (run_line)
}


