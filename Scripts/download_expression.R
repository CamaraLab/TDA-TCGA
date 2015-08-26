library(stringr)
library(getopt)
library(RCurl)
library(XML)

spec = matrix(c(
  "expression", "e",1, "character",
  # "project", "p", 1, "character",
  "url", "u",1,"character"
), byrow=TRUE, ncol=4)

arg<-getopt(spec) #Conmment this line for debug mode

url<-arg$url

html<-getURL(url)
html<-strsplit(html,"\n")
html<-unlist(html)
files<-html[grep (pattern = "*rsem.genes.results",html)]
files<-str_extract(files, ">unc.+\\.rsem.genes.results")
files<-substring(files,2)

dest_dir<-paste0(arg$expression)

if (Sys.info()['sysname']=="Windows") {
  method<-"auto"
  setInternet2(use = TRUE)} else {method<-"curl"}

sapply(files,function(file) {
  download.file(paste0(url,file),paste0(dest_dir,"/",file),method = method)
})
