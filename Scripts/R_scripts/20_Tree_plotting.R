#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


#-----------------------Figure 2 Plots------------------------------------------------------------------------------


#######################################################################
# Plot EF-Tu Tree with data that includes genera, similarity, and ROS induction
#######################################################################

EFTu_full_protein_tree <- read.tree("./../../Analyses/Protein_alignments_and_trees/EfTu/EFTu_full_length_alignment.treefile")
EFTu_full_protein_tree <- phangorn::midpoint(EFTu_full_protein_tree, node.labels='label')



for (i in 1:length(EFTu_full_protein_tree$tip.label)){
  EFTu_full_protein_tree$tip.label[i] <- paste(strsplit(EFTu_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                               strsplit(EFTu_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                               sep = "|")
}

# subset elf18 variant data
MAMP_elf18_hit_data <- subset(filtered_hold_MAMP_seqs, MAMP_Hit == "elf18_consensus")

# add info for ROS outcome 


MAMP_elf18_hit_data["New_label"] <- paste(MAMP_elf18_hit_data$Protein_Name, MAMP_elf18_hit_data$File_Name, sep = "|")
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[c("New_label", "Percent_Identity", "Genera")]
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[MAMP_elf18_hit_data$New_label %in% EFTu_full_protein_tree$tip.label,]


EFTu_tree <- ggtree(EFTu_full_protein_tree,  layout="rectangular", ladderize = T, size = 0.25, linetype = 1) %<+% MAMP_elf18_hit_data +
  geom_tippoint(aes(color = Genera), size = 0.05, alpha = 0.5, show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) +
  #layout_dendrogram() +
  #xlim(0, 0) 
  geom_treescale(x = -0.5, y = -3, linesize = 1, family = "Arial", offset = 28) +
  
  #ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.07, offset = 0.05, show.legend = FALSE) +
  
  #geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.1, offset = -1.06, show.legend = FALSE) +
  scale_fill_manual("Genera", values = Genera_colors) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.07, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  
  #geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.13,axis.params = list(line.color = "black"), offset = -2.17, show.legend = FALSE) +
  scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) 




# need to make data table to map bootstrap values onto tree
d <- EFTu_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 99,]

EFTu_tree <- EFTu_tree +
  geom_nodepoint(data=d,aes(label=label), shape = 21, size = 1, fill = "grey35") 










%>% ggtree::collapse(1396, 'mixed', fill="#7bc98f", color="#7bc98f")


elf18_tree + theme(legend.position = c(0.05, 0.3), legend.title = element_text("Genera", family = "Arial" ,
                                                                               face = "bold", color = "black", size = 12),
                   legend.text = element_text(family = "Arial", color = "black", size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))




# need to make data table to map bootstrap values onto tree
d <- elf18_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 70,]






##############################################
# load a tree file in newick tree format
##############################################


CSP_full_protein_tree <- read.tree("./../../Analyses/Protein_alignments_and_trees/cold_shock_protein/csp_domain_tree/CSD_hits_aligned.treefile")
CSP_full_protein_tree <- phangorn::midpoint(CSP_full_protein_tree, node.labels='label')


for (i in 1:length(CSP_full_protein_tree$tip.label)){
  CSP_full_protein_tree$tip.label[i] <- paste(strsplit(CSP_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                                strsplit(CSP_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                                sep = "|")
}


MAMP_csp22_hit_data <- subset(filtered_hold_MAMP_seqs, MAMP_Hit == "csp22_consensus")
MAMP_csp22_hit_data$MAMP_response <- NA

true_false_vector <- MAMP_csp22_hit_data$MAMP_Sequence %in% csp22_ROS_pos
MAMP_csp22_hit_data[which(true_false_vector == TRUE),9] <- "+"

true_false_vector <- MAMP_csp22_hit_data$MAMP_Sequence %in% csp22_ROS_neg
MAMP_csp22_hit_data[which(true_false_vector == TRUE),9] <- "-"

MAMP_csp22_hit_data["New_label"] <- paste(MAMP_csp22_hit_data$Protein_Name, MAMP_csp22_hit_data$File_Name, sep = "|")
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[c("New_label", "Percent_Identity", "Genera", "MAMP_response")]
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[MAMP_csp22_hit_data$New_label %in% CSP_full_protein_tree$tip.label,]



# equal_angle

CSP_tree <- ggtree(CSP_full_protein_tree, layout = "rectangular", ladderize = T, size = 0.18, linetype = 1) %<+% MAMP_csp22_hit_data +
  #geom_treescale(x = -0.3) +
  geom_treescale(x = -0.9, y = -3, linesize = 1, family = "Arial", offset = 45) +
  geom_tippoint(aes(color = Genera), size = 0.05, alpha = 0.5, show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) +
  #geom_tiplab(size = 0.1, align = T, linesize = 0.1) +
  #geom_tiplab(aes(label = Plant_Species), size = 1.8, color = "black", offset = 0.3, linesize = 0.1, align = T, family = "Arial", fontface = 'italic')
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.2, offset = 0.1, show.legend = FALSE) +
  #geom_tippoint(aes(color = Genera), size = 0.05, show.legend = FALSE) +
  scale_fill_manual("Genera", values = Genera_colors) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.2, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  
  #geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.13,axis.params = list(line.color = "black"), offset = -2.17, show.legend = FALSE) +
  scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) +
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = MAMP_response), width = 0.2, offset = 0.1) +
  scale_fill_manual("MAMP Response", values = c('white','black')) +
  theme(legend.direction = 'horizontal', legend.position = 'bottom')

CSP_tree <- CSP_tree +
  geom_tree(aes(color = Genera)) +
  scale_color_manual("Genera", values = Genera_colors) 
  



  #geom_nodelab(size = 7, col= "red") +
  #geom_tiplab(size = 0.3, align = T, linesize = 0.015)+
  #layout_dendrogram() +

  #theme(legend.position = "none") +
  #geom_treescale(fontsize = 2, linesize = 0.3, x = 1, y = 0) +

  #geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.07, offset = 0.05, show.legend = FALSE) +
  
  #scale_fill_manual("Genera", values = Genera_colors) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.2, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  
  #geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.13,axis.params = list(line.color = "black"), offset = -2.17, show.legend = FALSE) +
  scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) 





#csp22_tree <- csp22_tree + ggplot2::theme(legend.position = c(0.02, 0.3), legend.title = element_text("Genera", family = "Arial" ,
#                                                                               face = "bold", color = "black", size = 12),
#                   legend.text = element_text(family = "Arial", color = "black", size = 12)) +
#  guides(colour = guide_legend(override.aes = list(size = 3)))


# need to make data table to map bootstrap values onto tree
d <- csp22_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 70,]

csp22_tree <- csp22_tree + 
  geom_nodepoint(data = d, aes(label = label),size = 0.6, alpha = 0.6)


WP_csp22 <- strsplit(csp22_tree$data$label, "|", fixed = T)
for (i in 1:length(WP_csp22)) {
  WP_csp22[[i]] <- WP_csp22[[i]][1]
}


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



















##############################################
# subset 'core' gene phylogeny by genera
##############################################



Agro_tree %<+% map_data +
  geom_fruit(geom = geom_tile, mapping = aes(color = Genera, fill = Genera), width = 0.04,
             offset = -1.02, axis.params = list(line.color = "black")) +
  scale_color_manual("Genera", values = Genera_colors) +
  scale_fill_manual("Genera", values = Genera_colors) +
  theme(legend.position = "none") +
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
  scale_fill_gradient(low = "white", high = "black", breaks = c(0,2,4,6,8,10,12,14), 
                      guide = "legend")

#geom_hilight(data=nodedf, mapping=aes(node=node),
#             extendto=6.8, alpha=0.3, fill="grey", color="grey50",
#             size=0.05) +

#######################################################################
# plotting protein tree Genus Specific CSP tress
#######################################################################


# Clavibacter CSP
Clavibacter_CSP_protein_tree <- read.tree("./../Protein_alignments_and_trees/Genera_specific_trees/Clavibacter/CSPs/Clavibacter_CSPs_alignment.treefile")
Clavibacter_CSP_protein_tree <- phangorn::midpoint(Clavibacter_CSP_protein_tree, node.labels='label')
for (i in 1:length(Clavibacter_CSP_protein_tree$tip.label)){
  Clavibacter_CSP_protein_tree$tip.label[i] <- paste(strsplit(Clavibacter_CSP_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                                strsplit(Clavibacter_CSP_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                                sep = "|")
}
MAMP_csp22_hit_data <- subset(filtered_hold_MAMP_seqs, MAMP_Hit == "csp22_consensus")
CSP_gene_length <- All_target_by_annotation[All_target_by_annotation$Protein_Name %in% MAMP_csp22_hit_data$Protein_Name,]
CSP_gene_length <- CSP_gene_length[CSP_gene_length$Protein_Name %in% MAMP_csp22_hit_data$Protein_Name,]
MAMP_csp22_hit_data <- cbind(MAMP_csp22_hit_data, )
MAMP_csp22_hit_data["New_label"] <- paste(MAMP_csp22_hit_data$Protein_Name, MAMP_csp22_hit_data$File_Name, sep = "|")
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[c("New_label", "Percent_Identity", "Genera")]
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[MAMP_csp22_hit_data$New_label %in% Clavibacter_CSP_protein_tree$tip.label,]



Clavibacter_CSP_protein_tree_plot <- ggtree(Clavibacter_CSP_protein_tree, size = 0.2, layout = 'circular', ladderize = T) %<+% MAMP_csp22_hit_data +
  geom_treescale(fontsize = 2, linesize = 0.3, x = 2, y = 0) +
  geom_tiplab(size = 0.7, offset = 0.01) +
  geom_tippoint(aes(color = Percent_Identity), size = 1) +
  scale_color_viridis(option="magma", name = "Percent AA Similarity\n", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) +
  theme(legend.position = "none")
  

# Rhodococcus CSP
Rhodococcus_CSP_protein_tree <- read.tree("./../Protein_alignments_and_trees/Genera_specific_trees/Rhodococcus/CSPs/Rhodococcus_CSPs_alignment.treefile")
Rhodococcus_CSP_protein_tree <- phangorn::midpoint(Rhodococcus_CSP_protein_tree, node.labels='label')
for (i in 1:length(Rhodococcus_CSP_protein_tree$tip.label)){
  Rhodococcus_CSP_protein_tree$tip.label[i] <- paste(strsplit(Rhodococcus_CSP_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                                     strsplit(Rhodococcus_CSP_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                                     sep = "|")
}
MAMP_csp22_hit_data <- subset(filtered_hold_MAMP_seqs, MAMP_Hit == "csp22_consensus")
MAMP_csp22_hit_data["New_label"] <- paste(MAMP_csp22_hit_data$Protein_Name, MAMP_csp22_hit_data$File_Name, sep = "|")
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[c("New_label", "Percent_Identity", "Genera")]
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[MAMP_csp22_hit_data$New_label %in% Rhodococcus_CSP_protein_tree$tip.label,]



Rhodococcus_CSP_protein_tree_plot <- ggtree(Rhodococcus_CSP_protein_tree, size = 0.15, layout = 'circular', ladderize = T) %<+% MAMP_csp22_hit_data +
  geom_treescale(fontsize = 2, linesize = 0.3, x = -1, y = 0) +
  geom_tiplab(size = 0.4, offset = 0.05) +
  geom_tippoint(aes(color = Percent_Identity), size = 0.8) +
  scale_color_viridis(option="magma", name = "Percent AA Similarity\n", 
                      breaks = c(0, 20, 40, 60, 80, 100),
                      labels = c(0, 20, 40, 60, 80, 100),
                      limits = c(0,100)) +
  theme(legend.position = "none")

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
  theme_tree() +