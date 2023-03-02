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


csp_catagorization_tree <- function(tree_file, catagorization_file){
  tree_file <- phangorn::midpoint(tree_file, node.labels='label')
  
  for (i in 1:length(tree_file$tip.label)){
    tree_file$tip.label[i] <- paste(strsplit(tree_file$tip.label[i], "|", fixed = T)[[1]][1],
                                    strsplit(tree_file$tip.label[i], "|", fixed = T)[[1]][5],
                                    sep = "|")
  }
  
  catagorization_file["New_label"] <- paste(catagorization_file$WP_Tag, catagorization_file$File_name, sep = "|")
  catagorization_file <- catagorization_file[c("New_label", "Percent_Identity", "Type_Catagorization")]
  catagorization_file <- catagorization_file[catagorization_file$New_label %in% tree_file$tip.label,]
  
  
  save_tree <- ggtree(tree_file, size = 0.5, layout = 'rectangular', ladderize = T) %<+% catagorization_file +
    geom_treescale(fontsize = 4, linesize = 1.4, x = 0.2, y = 0, offset = 3) +
    geom_tippoint(aes(color = Type_Catagorization), size = 1) +
    
    ggnewscale::new_scale_fill() +
    geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.05, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
    scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                       breaks = c(0, 20, 40, 60, 80, 100),
                       labels = c(0, 20, 40, 60, 80, 100),
                       limits = c(0,100)) 
  
  return(save_tree)
}



# Clavibacter CSP
Clavibacter_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/Clavibacter/CSPs/reformat_clavi_cds_hits_aligned.tre")
import_clav_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Clavibacter")
csp_catagorization_tree(Clavibacter_CSP_protein_tree, import_clav_csp_types)


#  geom_treescale(fontsize = 4, linesize = 1.4, x = 0.2, y = 0, offset = 3) +





# Rhodococcus CSP
Rhodococcus_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/Rhodococcus/CSPs/reformat_rhodo_cds_hits_aligned.tre")
import_rhodo_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Rhodococcus")
csp_catagorization_tree(Rhodococcus_CSP_protein_tree, import_rhodo_csp_types)


#geom_treescale(fontsize = 4, linesize = 1.2, x = -1, y = 0, offset = 100) +
  



# Leifsonia CSP
Leifsonia_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/Leifsonia/reformat_leif_cds_hits_aligned.tre")
import_leif_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Leifsonia")
csp_catagorization_tree(Leifsonia_CSP_protein_tree, import_leif_csp_types)


#geom_treescale(fontsize = 4, linesize = 1.2, x = -1, y = 0, offset = 1) +




# Ralstonia CSP
Ralstonia_CSP_protein_tree <- read.tree(file = "./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/Ralstonia/reformat_rals_cds_hits_aligned.tre")
import_rals_csp_types <- xlsx:::read.xlsx("./../../Analyses/Catagorizing_MAMPs/Genera_specific_analysis/MAMP_Catagorization.xlsx", sheetName = "Ralstonia")




filepath <- system.file("examples", "meme.xml", package = "ggmotif")
motif_extract <- getMotifFromMEME(data = file.choose(), format="xml")
motif_plot <- motifLocation(data = motif_extract)
motif_plot +
  ggsci::scale_fill_aaas() +
  theme(axis.text.y = element_text(size = 3))




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







