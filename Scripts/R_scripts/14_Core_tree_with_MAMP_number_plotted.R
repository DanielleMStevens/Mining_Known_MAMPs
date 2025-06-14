#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/08/2023
# Script Purpose: Plotting MAMP hits onto phylogenetic tree
# Inputs: GToTree tre file output
# Outputs: phylogenetic tree 
#-----------------------------------------------------------------------------------------------


##############################################
# phylogenetic tree based on 74 housekeeping genes
##############################################

# import number of copies of each mamp on a per genome basis
map_data <- datasettable[,c(6,1:5,7)]

# filter genomes based on ANI analysis
map_data <- map_data[!map_data$Filename %in% Genomes_to_check$Genome1,]


# Note: hold_copy_number is from final list in script 12_Finalize_MAMP_list_post_ANI.R
csp22_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "csp22_consensus")
map_data <- cbind(map_data, "csp22_consensus" = csp22_copy_number[match(map_data[,1], csp22_copy_number[,2]),5])

elf18_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "elf18_consensus")
map_data <- cbind(map_data, "elf18_consensus" = elf18_copy_number[match(map_data[,1], elf18_copy_number[,2]),5])
map_data[is.na(map_data$elf18_consensus),9] <- 0

flg22_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "flg22_consensus")
map_data <- cbind(map_data, "flg22_consensus" = flg22_copy_number[match(map_data[,1], flg22_copy_number[,2]),5])
map_data[is.na(map_data$flg22_consensus),10] <- 0

flg28_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "flgII-28")
map_data <- cbind(map_data, "flgII-28" = flg28_copy_number[match(map_data[,1], flg28_copy_number[,2]),5])
map_data[is.na(map_data$`flgII-28`),11] <- 0

nlp20_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "nlp20_consensus")
map_data <- cbind(map_data, "nlp20_consensus" = nlp20_copy_number[match(map_data[,1], nlp20_copy_number[,2]),5])
map_data[is.na(map_data$`nlp20_consensus`),12] <- 0

#map_data_to_remove <- map_data[!map_data$Filename %in% filtered_hold_MAMP_seqs$File_Name,]
#map_data <- map_data[!map_data$Filename %in% map_data_to_remove$Filename,]

core_gene_phylo <- read.tree("./../../Analyses/All_bacteria_phylogenomic_tree/Figure1B_Phylogenetic_Tree/GToTree_output/GToTree_output.tre")
core_gene_phylo <- phangorn::midpoint(core_gene_phylo, node.labels='label')

for (i in 1:length(core_gene_phylo$tip.label)){
  accession_label <- paste(strsplit(core_gene_phylo$tip.label[[i]], "_")[[1]][1],
                           strsplit(core_gene_phylo$tip.label[[i]], "_")[[1]][2],
                           sep = "_")
  core_gene_phylo$tip.label[[i]] <- datasettable[datasettable$Assembly_Accession %in% accession_label,6]
}


# plot full tree -> basic
full_tree <- ggtree::ggtree(core_gene_phylo,  ladderize = T, size = 0.09, linetype = 1) 


# add details to main tree regarding genome generas and MAMP abundance 
full_tree <- full_tree %<+% map_data +
  geom_fruit(geom = geom_tile, mapping = aes(color = Genera, fill = Genera), width = 0.04,
             offset = -1.02, axis.params = list(line.color = "black"), show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) +
  scale_fill_manual("Genera", values = Genera_colors) +
  layout_dendrogram() +
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(fill = csp22_consensus), width = 0.04,
             offset = -1.08) +
  geom_fruit(geom = geom_tile, mapping = aes(fill = elf18_consensus), width = 0.04, 
             offset = -2.2) +
  geom_fruit(geom = geom_tile, mapping = aes(fill = flg22_consensus), width = 0.04,
             offset = -4.44) + 
  geom_fruit(geom = geom_tile, mapping = aes(fill = `flgII-28`), width = 0.04,
             offset = -8.92) +
  geom_fruit(geom = geom_tile, mapping = aes(fill = nlp20_consensus), width = 0.04,
             offset = -17.88) +
  
  # ignore - other ways to gradient fill in mamp abundance 
  #scale_fill_stepsn(colours = c("#ffffff","#e5e5e5","#c7c7c7","#999999","black"),
  #                  values = scales::rescale(c(0,1,2,4,8), from = c(0,8)),
  #                  breaks = c(0,1,2,4,8),
  #                limits = c(0,8),
  #                guide = "legend", name = "MAMP Abundance") +
  
  #scale_fill_stepsn(
  #  colours=c("#ffffff","#e5e5e5","#c7c7c7","#999999","black"),
  #  values=c(  0,       1,  2,    4,      8)/8) +

  scale_fill_gradient(low = "white", high = "black", breaks = c(0,1,2,4,8), limits = c(0,15), 
                      guide = "legend", name = "MAMP Abundance") +
  
  theme(legend.direction = "horizontal", 
        legend.position = "bottom")



#ggsave(full_tree, filename = "./../../Figures/Figure_1/Full_phylogenomic_tree_with_MAMPs_v2.pdf", device = cairo_pdf, width = 7, height = 3.2, units = "in")


#---------------------------------Notes for Figure 1--------------------------------
# The PDF of the plotted tree was exported and imported (via cairo pockage) in insckape
# The following asthetic modificaitons were made: added black outline around each mamp
# abundance plot. Also added an outline around the bar for genera coloring. Finally added
# a 




#### ---------------------------------for MPMI Presentation---------------------------------
#full_tree <- ggtree(core_gene_phylo, layout="daylight", branch.length = 'none',  ladderize = T, size = 0.12, linetype = 1)
#full_tree <- ggtree::ggtree(core_gene_phylo, layout = 'fan', open.angle = 110, ladderize = T, size = 0.15, linetype = 1) 

#full_tree %<+% map_data +
#  geom_fruit(geom = geom_tile, mapping = aes(color = Genera, fill = Genera), width = 0.04,
#             axis.params = list(line.color = "black"), show.legend = FALSE) +
#  scale_color_manual("Genera", values = Genera_colors) +
  scale_fill_manual("Genera", values = Genera_colors) 


