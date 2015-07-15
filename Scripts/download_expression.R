url<-"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/gbm/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/unc.edu_GBM.IlluminaHiSeq_RNASeqV2.Level_3.1.2.0/"
html<-readLines(url)
files<-html[grep (pattern = "*rsem.genes.results",code)]
files<-str_extract(files, ">unc.+\\.rsem.genes.results")
files<-substring(files,2)
files[1]
sapply(files,function(file) {
  download.file(paste0(url,file),file)
})
