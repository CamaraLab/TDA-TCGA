library(maftools)


# CREATING INPUT MUTATION TABLE -------------------------------------------

setwd("/home/rstudio/documents/Messy_test_data/LGG_maf_files")
files <- list.files(pattern = "TCGA-*")

syn_muts <- data.frame("Sample"= NA,"Gene" = NA,"Mutation" = NA, "Type" = NA)
nonsyn_muts <- syn_muts

for(maf in files){
  update <- read.maf(read.delim(maf,comment.char="#"))
  nonsyn <- update@data
  Sample <- substr(nonsyn$Tumor_Sample_Barcode, 1,16)
  Gene <- paste0(nonsyn$Hugo_Symbol,".",nonsyn$Entrez_Gene_Id)
  Mutation <- nonsyn$Protein_Change
  Type <- nonsyn$Variant_Classification
  nonsyn_muts <- rbind(nonsyn_muts,data.frame(Sample,Gene,Mutation,Type))
  
  
  syn <- update@maf.silent
  Sample <- substr(syn$Tumor_Sample_Barcode, 1,16)
  Gene <- paste0(syn$Hugo_Symbol,".",syn$Entrez_Gene_Id)
  Mutation <- syn$Protein_Change
  Type <- syn$Variant_Classification
  syn_muts <- rbind(syn_muts,data.frame(Sample,Gene,Mutation,Type))

}

# Get rid of first 'NA' row
nonsyn_muts <- nonsyn_muts[2:nrow(nonsyn_muts),]
row.names(nonsyn_muts) <- seq(1,nrow(nonsyn_muts),1)

syn_muts <- syn_muts[2:nrow(syn_muts),]
row.names(syn_muts) <- seq(1,nrow(syn_muts),1)

# Full Mutation List
LGG_muts <- rbind(nonsyn_muts,syn_muts)

write.csv(nonsyn_muts,"LGG_NonSynMuts.txt")
write.csv(syn_muts,"LGG_SynMuts.txt")
write.csv(LGG_muts,"LGG_Muts.txt")
