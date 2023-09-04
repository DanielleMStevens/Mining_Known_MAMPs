#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


#######################################################################
# plotting protein tree Genus Specific CSP tress
#######################################################################


csp_catagorization_tree <- function(tree_file, catagorization_file, tree_collapse_adj, ring_thickness, offset_value){
  tree_file <- phangorn::midpoint(tree_file, node.labels='label')
  
  for (i in 1:length(tree_file$tip.label)){
    tree_file$tip.label[i] <- paste(strsplit(tree_file$tip.label[i], "|", fixed = T)[[1]][1],
                                    strsplit(tree_file$tip.label[i], "|", fixed = T)[[1]][5],
                                    sep = "|")
  }
  
  catagorization_file["New_label"] <- paste(catagorization_file$WP_Tag, catagorization_file$File_name, sep = "|")
  catagorization_file <- catagorization_file[c("New_label", "Percent_Identity", "Type_Catagorization")]
  catagorization_file <- catagorization_file[catagorization_file$New_label %in% tree_file$tip.label,]
  
  
  save_tree <- ggtree(tree_file, size = 0.5, layout = 'circular', ladderize = T) %<+% catagorization_file +
    geom_treescale(fontsize = 3.5, linesize = 1, x = tree_collapse_adj, y = 0, offset = offset_value) +
    geom_tippoint(aes(color = Type_Catagorization), size = 1, show.legend = TRUE) +
    
    ggnewscale::new_scale_fill() +
    geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = ring_thickness, offset = 0.2, axis.params = list(line.color = "black"), show.legend = FALSE) +
    scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                       breaks = c(0, 20, 40, 60, 80, 100),
                       labels = c(0, 20, 40, 60, 80, 100),
                       limits = c(0,100)) 

  return(save_tree)
}



# -----------------------------Clavibacter CSP--------------------------------------------------------
Clavibacter_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_CSPs/CSP_types/Clavibacter/reformat_clavi_domain.tre")
import_clav_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Clavibacter")
motif_extract <- getMotifFromMEME(data = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Clavibacter/Clavi_meme_results.xml", format="xml")


motif_test <- motif_extract
motif_test$input.seq.id <- as.character(motif_test$input.seq.id)
for (i in 1:nrow(motif_test)){
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "csp22_consensus\\|Full_Seq\\|Clavibacter\\|","")
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "\\|[0-9]","")
}
motif_test <- motifLocation(motif_test[!duplicated(motif_test[c(2,3,5,6,7,8,9)]),])

# csp meme results
motif_test +
  ggsci::scale_fill_aaas() +
  ylab("") +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"))


csp_catagorization_tree(Clavibacter_CSP_protein_tree, import_clav_csp_types, -0.6, 0.11, 0.2)





# -----------------------------Rhodococcus CSP--------------------------------------------------------
Rhodococcus_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Rhodococcus/CSPs/reformat_rhodo_cds_hits_aligned.tre")

Rhodococcus_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_CSPs/CSP_types/Rhodococcus/reformat_rhoco_domain.tre")


import_rhodo_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Rhodococcus")
motif_extract <- getMotifFromMEME(data = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Rhodococcus/Rhodo_meme_results.xml", format="xml")


motif_test <- motif_extract
motif_test$input.seq.id <- as.character(motif_test$input.seq.id)
for (i in 1:nrow(motif_test)){
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "csp22_consensus\\|Full_Seq\\|Rhodococcus\\|","")
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "\\|[0-9]","")
}
motif_test <- motifLocation(motif_test[!duplicated(motif_test[c(2,3,5,6,7,8,9)]),])

# csp tree
csp_catagorization_tree(Rhodococcus_CSP_protein_tree, import_rhodo_csp_types, -0.8, 0.3, 100)

# csp meme results
motif_test +
  ggsci::scale_fill_aaas() +
  ylab("") +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"))




# -----------------------------Leifsonia CSP--------------------------------------------------------
Leifsonia_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Leifsonia/CSPs/reformat_leif_cds_hits_aligned.tre")
import_leif_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Leifsonia")
csp_catagorization_tree(Leifsonia_CSP_protein_tree, import_leif_csp_types, -0.8, 0.2, 10)

motif_extract <- getMotifFromMEME(data = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Leifsonia/Leif_meme_results.xml", format="xml")


motif_test <- motif_extract
motif_test$input.seq.id <- as.character(motif_test$input.seq.id)
for (i in 1:nrow(motif_test)){
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "csp22_consensus\\|Full_Seq\\|Leifsonia\\|","")
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "\\|[0-9]","")
}
motif_test <- motifLocation(motif_test[!duplicated(motif_test[c(2,3,5,6,7,8,9)]),])

# csp meme results
motif_test +
  ggsci::scale_fill_aaas() +
  ylab("") +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"))




# -----------------------------Erwinia CSP--------------------------------------------------------
Erwinia_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Erwinia/reformat_erwin_cds_hits_aligned.tre")
import_erwin_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Erwinia")
csp_catagorization_tree(Erwinia_CSP_protein_tree, import_erwin_csp_types, -0.8, 0.2, 10)

motif_extract <- getMotifFromMEME(data = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Erwinia/Erwin_meme_results.xml", format="xml")


motif_test <- motif_extract
motif_test$input.seq.id <- as.character(motif_test$input.seq.id)
for (i in 1:nrow(motif_test)){
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "csp22_consensus\\|Full_Seq\\|Erwinia\\|","")
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "\\|[0-9]","")
}
motif_test <- motifLocation(motif_test[!duplicated(motif_test[c(2,3,5,6,7,8,9)]),])

# csp meme results
motif_test +
  ggsci::scale_fill_aaas() +
  ylab("") +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"))




# -----------------------------Xanthomonas CSP--------------------------------------------------------
Xanthomonas_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Xanthomonas/reformat_xanth_domain.tre")
import_xanth_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Xanthomonas")
csp_catagorization_tree(Xanthomonas_CSP_protein_tree, import_xanth_csp_types, -0.4, 0.2, 10)


motif_extract <- getMotifFromMEME(data = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Xanthomonas/Xanth_meme_results.xml", format="xml")


motif_test <- motif_extract
motif_test$input.seq.id <- as.character(motif_test$input.seq.id)
for (i in 1:nrow(motif_test)){
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "csp22_consensus\\|Full_Seq\\|Xanthomonas\\|","")
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "\\|[0-9]","")
}
motif_test <- motifLocation(motif_test[!duplicated(motif_test[c(2,3,5,6,7,8,9)]),])

# csp meme results
motif_test +
  ggsci::scale_fill_aaas() +
  ylab("") +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"))





# -----------------------------Ralstonia CSP--------------------------------------------------------
Ralstonia_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Ralstonia/reformat_rals_cds_hits_aligned.tre")
import_rals_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Ralstonia")
csp_catagorization_tree(Ralstonia_CSP_protein_tree, import_rals_csp_types, -0.8, 0.2, 10)

motif_extract <- getMotifFromMEME(data = "./../../Analyses/Catagorizing_CSPs/Genera_specific_analysis/Leifsonia/Leif_meme_results.xml", format="xml")


motif_test <- motif_extract
motif_test$input.seq.id <- as.character(motif_test$input.seq.id)
for (i in 1:nrow(motif_test)){
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "csp22_consensus\\|Full_Seq\\|Erwinia\\|","")
  motif_test$input.seq.id[i] <- str_replace(motif_test$input.seq.id[i], "\\|[0-9]","")
}
motif_test <- motifLocation(motif_test[!duplicated(motif_test[c(2,3,5,6,7,8,9)]),])

# csp meme results
motif_test +
  ggsci::scale_fill_aaas() +
  ylab("") +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"))
















 #######################################################################
# plotting tree for elf18
#######################################################################


elf18_tree <- read.tree("./../Protein_alignments_and_trees/EfTu/MAMP_fasta/elf18_alignment.treefile")
elf18_tree <- phangorn::midpoint(elf18_tree, node.labels='label')



for (i in 1:length(elf18_tree$tip.label)){
  elf18_tree$tip.label[i] <- paste(strsplit(elf18_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                   strsplit(elf18_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                   sep = "|")
}


MAMP_elf18_hit_data <- subset(hold_MAMP_seqs, MAMP_Hit == "elf18_consensus")
MAMP_elf18_hit_data["New_label"] <- paste(MAMP_elf18_hit_data$Protein_Name, MAMP_elf18_hit_data$File_Name, sep = "|")
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[c("New_label", "Percent_Identity", "Genera")]
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[MAMP_elf18_hit_data$New_label %in% elf18_tree$tip.label,]


ggtree(elf18_tree, layout="circular", ladderize = T, size = 0.35, linetype = 1) %<+% MAMP_elf18_hit_data +
  geom_tippoint(aes(color = Genera), size = 0.2) +
  scale_color_manual("Genera", values = Genera_colors) +
  theme(legend.position = "none") +
  geom_treescale(fontsize = 2, linesize = 0.3, x = 1.5, y = 0) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.1) +
  scale_fill_manual("Genera", values = Genera_colors) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.2,
             offset = 0.08, axis.params = list(line.color = "black")) +
  scale_fill_viridis(option="magma", name = "Percent AA Similarity\n", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) 




#######################################################################
# plotting tree for flg22
#######################################################################


flg22_tree <- read.tree("./../Protein_alignments_and_trees/Flagellin/MAMP_fasta/flg22_alignment.treefile")
flg22_tree <- phangorn::midpoint(flg22_tree, node.labels='label')



for (i in 1:length(flg22_tree$tip.label)){
  flg22_tree$tip.label[i] <- paste(strsplit(flg22_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                   strsplit(flg22_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                   sep = "|")
}


MAMP_flg22_hit_data <- subset(hold_MAMP_seqs, MAMP_Hit == "flg22_consensus")
MAMP_flg22_hit_data["New_label"] <- paste(MAMP_flg22_hit_data$Protein_Name, MAMP_flg22_hit_data$File_Name, sep = "|")
MAMP_flg22_hit_data <- MAMP_flg22_hit_data[c("New_label", "Percent_Identity", "Genera")]
MAMP_flg22_hit_data <- MAMP_flg22_hit_data[MAMP_flg22_hit_data$New_label %in% flg22_tree$tip.label,]


ggtree(flg22_tree, layout="circular", ladderize = T, size = 0.35, linetype = 1) %<+% MAMP_flg22_hit_data +
  geom_tippoint(aes(color = Genera), size = 0.5) +
  scale_color_manual("Genera", values = Genera_colors) +
  theme(legend.position = "none") +
  geom_treescale(fontsize = 2, linesize = 0.3, x = 1.5, y = 0) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.1) +
  scale_fill_manual("Genera", values = Genera_colors) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.2,
             offset = 0.08, axis.params = list(line.color = "black")) +
  scale_fill_viridis(option="magma", name = "Percent AA Similarity\n", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) 


############# for rectangular treee
#max_value_edge <- max(elf18_full_protein_tree$edge.length)*2.5
theme_tree2() + #this next three lines are a hacky way to "rescale" the tree 
  xlim(0, max_value_edge)+
  theme_tree() 






#######################################################################
# plotting protein tree for Flagellin
#######################################################################


flg22_full_protein_tree <- read.tree("./../Protein_alignments_and_trees/Flagellin/filC_full_length_alignment.treefile")
flg22_full_protein_tree <- phangorn::midpoint(flg22_full_protein_tree, node.labels='label')



for (i in 1:length(flg22_full_protein_tree$tip.label)){
  flg22_full_protein_tree$tip.label[i] <- paste(strsplit(flg22_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                                strsplit(flg22_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                                sep = "|")
}


MAMP_flg22_hit_data <- subset(hold_MAMP_seqs, MAMP_Hit == "flg22_consensus")
MAMP_flg22_hit_data["New_label"] <- paste(MAMP_flg22_hit_data$Protein_Name, MAMP_flg22_hit_data$File_Name, sep = "|")
MAMP_flg22_hit_data <- MAMP_flg22_hit_data[c("New_label", "Percent_Identity", "Genera")]
MAMP_flg22_hit_data <- MAMP_flg22_hit_data[MAMP_flg22_hit_data$New_label %in% flg22_full_protein_tree$tip.label,]


ggtree(flg22_full_protein_tree, layout = "circular", ladderize = T, size = 0.35, linetype = 1) %<+% MAMP_flg22_hit_data +
  geom_tippoint(aes(color = Genera), size = 0.8, show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) +
  #layout_dendrogram() +
  
  #theme(legend.position = "none") +
  #geom_treescale(fontsize = 2, linesize = 0.3, x = 2, y = 0) +
  ggnewscale::new_scale_fill() +
  #geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.1, offset = -1.06, show.legend = FALSE) +
  #scale_fill_manual("Genera", values = Genera_colors) +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.2, offset = 0.04, show.legend = FALSE) +
  scale_fill_manual("Genera", values = Genera_colors) +
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.2, offset = 0.04, axis.params = list(line.color = "black"), show.legend = FALSE) +
  
  #geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.13, axis.params = list(line.color = "black"), offset = -2.17, show.legend = FALSE) +
  scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) 



flg22_tree <- flg22_tree + theme(legend.position = c(0.05, 0.3), legend.title = element_text("Genera", family = "Arial" ,
                                                                                             face = "bold", color = "black", size = 12),
                                 legend.text = element_text(family = "Arial", color = "black", size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))


ggsave(flg22_tree, 
       filename = "./Figures/Plot_Flagellin_protein_tree.pdf",
       device = cairo_pdf, width = 12, height = 9, units = "in")






#----------------------------------------------------------------------------------------------------------
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


csp_protein_clusters <- data.frame(read_tsv(file = "./../../Analyses/Protein_alignments_and_trees/cold_shock_protein/csp_clustering/clusterRes_cluster.tsv"))
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












