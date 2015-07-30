PROJECT_NAME<-"LUSC"
wd<-paste0("c:/Users/Udi/Documents/TCGA-DATA/",PROJECT_NAME)
setwd(wd)
url<-"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/lusc/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.Level_3.1.9.0/"
html<-readLines(url)
files<-html[grep (pattern = "*rsem.genes.results",html)]
files<-str_extract(files, ">unc.+\\.rsem.genes.results")
files<-substring(files,2)
sapply(files,function(file) {
  download.file(paste0(url,file),paste0("./Expression/",file))
})
