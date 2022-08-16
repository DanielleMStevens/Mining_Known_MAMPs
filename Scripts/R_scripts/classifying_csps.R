#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------

######################################################################
# load output from MMseqs to classify protein clusters of csp proteins
######################################################################


csp_protein_clusters <- data.frame(read_tsv(file = "./../Protein_alignments_and_trees/cold_shock_protein/csp_clustering/clusterRes_cluster.tsv"))
colnames(csp_protein_clusters) <- c("cluster protein reference", "all proteins within cluster")


hold_MAMP_score <- list()
hold_genera_ID <- list()
for (i in 1:nrow(csp_protein_clusters)){
  pull_protein_tag <- str_split(csp_protein_clusters$`all proteins within cluster`[i], "\\|")[[1]][1]
  match_by_tag <- hold_MAMP_seqs[hold_MAMP_seqs$Protein_Name %in% pull_protein_tag,]
  pull_file_name <- str_split(csp_protein_clusters$`all proteins within cluster`[i], "\\|")[[1]][5]
  
  # pull score
  if(nrow(match_by_tag[match_by_tag$File_Name %in% pull_file_name,]) == 1){
    hold_value <- match_by_tag[match_by_tag$File_Name %in% pull_file_name,3]
    hold_genera <- match_by_tag[match_by_tag$File_Name %in% pull_file_name,7]
  }
  
  # pull score of first iteractions - need to find why therea re dublets
  if(nrow(match_by_tag[match_by_tag$File_Name %in% pull_file_name,]) > 1){
    hold_value <- match_by_tag[match_by_tag$File_Name %in% pull_file_name,]
    hold_genera <- hold_value[1,7]
    hold_value <- hold_value[1,3]
  }
  hold_genera_ID[[i]] <- hold_genera
  hold_MAMP_score[[i]] <- hold_value
}


csp_protein_clusters <- cbind(csp_protein_clusters, "AA Score" = unlist(hold_MAMP_score))
csp_protein_clusters <- cbind(csp_protein_clusters, "Genera" = unlist(hold_genera_ID))


n_number <- as.data.frame(csp_protein_clusters %>% group_by(`cluster protein reference`) %>% summarise(n=n()))
filter_csp_clusters <- subset(n_number, n_number$n > 9)

ggplot(csp_protein_clusters[csp_protein_clusters$`cluster protein reference` %in% filter_csp_clusters$`cluster protein reference`,], 
       aes(x = `cluster protein reference`, y = `AA Score`, color = Genera)) +
  geom_violinhalf(color = "black", trim = T, scale = "width") +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.15) +
  geom_jitter(height = 1.5, size =0.2)+
  scale_color_manual("Genera", values = Genera_colors) +
  my_ggplot_theme +
  xlab("Clusters\n") +
  ylab("\nPercent AA Similarity") +
  scale_y_continuous(limits = c(0,100),
                     breaks = c(0,20,40,60,80,100)) +
  #ggtitle(name_of_MAMP) +
  theme(legend.position = "none",panel.grid.major.x = element_line(size = 0.5, color = "grey92"))+
  coord_flip()



# this takes a long time to complete -1.5-2 hours
#pb <- txtProgressBar(min = 0, max = nrow(csp_full_length), style = 3)
#csp_sim_score <- data.frame("Protein_seq_1" = character(0), "Protein_seq_2" = character(0), "Percent_sim" = numeric(0))
#for (i in 1:nrow(csp_full_length)){
#  for (j in 2:nrow(csp_full_length)){
#    test2 <- Biostrings::pairwiseAlignment(csp_full_length$Sequence[[i]], csp_full_length$Sequence[[j]], type = "global", substitutionMatrix = BLOSUM62)
#    hold_score <- Biostrings::pid(test2, type = "PID1")
#    csp_sim_score <- rbind(csp_sim_score, data.frame("Protein_seq_1" = csp_full_length$Locus_Tag_Name[[i]],
#                                                     "Protein_seq_2" = csp_full_length$Locus_Tag_Name[[j]],
#                                                      "Percent_sim" = hold_score))
#  }
#  setTxtProgressBar(pb, i)
#  print(i)
#} 
