library(RayleighSelection)

compute_mut_load <- function(mut_table,nonsyn_muts=NA,syn_muts=NA,num_cores=1) {

  samples <- sort(unique(mut_table$Sample))
  gene_names <- sort(unique(mut_table$Gene))
  
  if(is.na(nonsyn_muts) || is.na(syn_muts)){
    nonsyn_type <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site",      
                     "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins",
                     "In_Frame_Ins", "Nonstop_Mutation")
    
    nonsyn_muts <- mut_table[mut_table$Type %in% nonsyn_type,]
    syn_muts <- mut_table[!(mut_table$Type %in% nonsyn_type),]
  }
  
  nonsyn_bin <- matrix(0,length(samples),length(gene_names)) %>% as.data.frame
  dimnames(nonsyn_bin) <- list(samples,gene_names)
  syn_bin <- nonsyn_bin
  
  t_nonsyn <- with(nonsyn_muts,table(Sample,Gene))
  nonsyn_bin[rownames(t_nonsyn),colnames(t_nonsyn)] <- t_nonsyn
  t_syn <- with(syn_muts,table(Sample,Gene,Type))
  syn_bin[rownames(t_syn),colnames(t_syn)] <- t_syn
  
  # nonsyn_bin<-ifelse(nonsyn_bin>0,1,0)
  # syn_bin <- ifelse(syn_bin>0,1,0)
  
  mutload <- rowSums(syn_bin+nonsyn_bin)
  hist(log10(mutload),breaks=100,plot=TRUE,main="LGG Mutational Load")


  missing_genes <- c("TCGA-DU-7014-01A","TCGA-DU-A7TI-01A","TCGA-HW-7493-01A","TCGA-TQ-A7RK-02B")
  newmat <- matrix(0,1,4)
  colnames(newmat) <- missing_genes
  mutload_mat <- cbind(t(mutload) %>% as.data.frame, newmat)
  mutload_mat <- mutload_mat[,order(colnames(mutload_mat))]
  pvals <- rayleigh_selection(gg,mutload_mat,num_cores=10)

}
