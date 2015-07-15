#setwd(".")
#setwd("C:/Users/Udi/Downloads/LUAD_3.1.14.0")
r=getOption("repos")
r["CRAN"]="http://cran.uk.r-project.org"
install.packages('igraph',repos=r)
install.packages('rgexf',repos=r)
install.packages('jsonlite',repos=r)
install.packages('getopt',repos=r)
install.packages("rhdf5",repos=r)
install.packages("data.table",repos=r)

