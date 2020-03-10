library(TDAmapper)
library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)
library(umap)
library(bioDist)

# INPUT AND CLEAN ---------------------------------------------------------

exp_table <- readline("Input gene expression table: ")
#exp_table <- (read.csv("C:/Users/Adam/Desktop/Camara Lab/TDA_TCGA_Project/Test Data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))
mut_table <- readline("Input mutation table: ")

#exp_table = t(exp_table)[1:1000,]

exp_table <- read.csv(exp_table,row.names=1,header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?"))) #row and column names from file
mut_table <- read.csv(mut_table,row.names=1,header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))

'''
colnames(df3)[1] <- c("gene")
colnames(df3)[2:ncol(df3)] <- paste("cell", 1:(ncol(df3)-1), sep="_")
colnames(mut_table) <- c("gene", "mut", "mut_type")
'''

exp_table <- exp_table[!duplicated(rownames(exp_table)),!duplicated(colnames(exp_table))] #No duplicated genes/samples
#exp_table <- exp_table[!duplicated(as.list(exp_table))]
mut_table <- mut_table[!duplicated(rownames(mut_table)),!duplicated(colnames(mut_table))]

if(all(exp_table$gene %in% mut_table$gene)) { #Gene names are same CHANGE
  if(all(mut_table$gene %in% exp_table$gene)) {
    print("Creating Distance Matrix...")
  }
} else {
  print("Unmatched genes between expression and mutation tables")
}


# MAPPER ------------------------------------------------------------------

plot.umap = function(emb) {
  ggplot(emb) +
    geom_point(aes(x=V1, y=V2, color=0), size=0.9)
}

build.mapper = function(dist, umap) {
  mapper2D(dist, umap, c(20, 20), percent_overlap=50, num_bins_when_clustering=10)
}

dist_matrix <- as.matrix(cor.dist(as.matrix(exp_table))) #from bioDist library in Bioconductor

umap_emb = umap(dist_matrix)$layout %>% as.data.frame
plot.umap(umap_emb)

mapperObj <- build.mapper(dist_matrix,umap_emb)
adj <- graph.adjacency(mapperObj$adjacency, mode="undirected")
plot(adj, layout = layout.auto(adj))


#### TRY OUT


