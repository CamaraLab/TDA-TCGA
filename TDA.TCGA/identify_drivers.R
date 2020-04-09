identify_drivers <- function(exp_table, mut_object, mutload_object, freq_threshold = 0.02, top_nonsyn_fraction = 350, permutations = 5000, q_threshold = 0.15,) {
  
  exp_table <- read.csv("/home/rstudio/documents/Messy_test_data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F)
  
  samples <- rownames(exp_table)

  # Freq threshold
  
  # NonSynonymous fraction threshold
  
  # Assess localization of genes passing threshold
  
  # Keep genes with p < 0.05 and q < 0.15
  
  
  
  
  
  
  
}