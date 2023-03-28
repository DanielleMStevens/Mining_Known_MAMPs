# kmer hasing

#sourmash sketch protein -p ksize 9 csp_full_length.fasta --singleton
#sourmash compare *.sig -o distances.cmp -k 9 --protein -p 10 --csv compare_csp_hash.csv

# import cvs file
kmer_hash <- read.csv(file.choose())
rownames(kmer_hash) <- colnames(kmer_hash)
kmer_hash_matrix <- as.matrix(kmer_hash)

# restructure dataframe such that comaprisons are based on 3 column format. 
# Jacard index less thna 0.1 will be filtered out and then exported before trying to load into cytoscape

kmer_hash_restructured <- data.frame("node" = character(0), "edge" = character(0), "jacard_index" = numeric(0))

for (i in 1:length(colnames(kmer_hash_matrix))){
  for (j in 1:length(rownames(kmer_hash_matrix))){
    if(kmer_hash_matrix[j,i] > 0){
      kmer_hash_restructured <- rbind(kmer_hash_restructured, data.frame("node" = colnames(kmer_hash_matrix)[i], 
                                                                         "edge" = rownames(kmer_hash_matrix)[j], 
                                                                         "jacard_index" = kmer_hash_matrix[j,i]))
      
    }
  }
}
  





#colnames(kmer_hash_matrix) <- colnames(kmer_hash)


ComplexHeatmap::Heatmap(matrix = kmer_hash_matrix,
                        
                        row_dend_reorder = T,
                        column_dend_reorder = T,
                        
                        show_column_names = F, 
                        show_row_names = F)


