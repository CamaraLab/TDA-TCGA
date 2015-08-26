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

#PROJECT_NAME<-arg$project
url<-arg$url
#url<-"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/prad/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/unc.edu_PRAD.IlluminaHiSeq_RNASeqV2.Level_3.1.14.0/"
#wd<-paste0(arg$folder,"/",PROJECT_NAME)
#setwd(wd)

#setInternet2(use = TRUE)
#url<-"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/skcm/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/unc.edu_SKCM.IlluminaHiSeq_RNASeqV2.Level_3.1.14.0/"
#PROJECT_NAME<-"SKCM"
#html<-readLines(url)
html<-getURL(url)
html<-strsplit(html,"\n")
html<-unlist(html)
files<-html[grep (pattern = "*rsem.genes.results",html)]
files<-str_extract(files, ">unc.+\\.rsem.genes.results")
files<-substring(files,2)

dest_dir<-paste0(arg$expression)
#dest_dir<-getwd()

if (Sys.info()['sysname']=="Windows") {method<-"auto"} else {method<-"curl"}

sapply(files,function(file) {
  download.file(paste0(url,file),paste0(dest_dir,"/",file),method = method)
})
