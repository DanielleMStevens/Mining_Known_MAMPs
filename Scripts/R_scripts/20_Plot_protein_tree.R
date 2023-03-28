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
MAMP_elf18_hit_data["Immunogenicity"] <- "NT"
for (i in 1:nrow(elf18_ros_screen_outcome)){
  MAMP_elf18_hit_data[MAMP_elf18_hit_data$MAMP_Sequence %in% elf18_ros_screen_outcome$Sequence[i],9] <- elf18_ros_screen_outcome$Outcome[i]
} 


# add in new label to match with protein tree tips
MAMP_elf18_hit_data["New_label"] <- paste(MAMP_elf18_hit_data$Protein_Name, MAMP_elf18_hit_data$File_Name, sep = "|")
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[c("New_label", "Percent_Identity", "Genera", "Immunogenicity")]
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[MAMP_elf18_hit_data$New_label %in% EFTu_full_protein_tree$tip.label,]


EFTu_tree <- ggtree(EFTu_full_protein_tree,  layout="rectangular", ladderize = T, size = 0.3, linetype = 1) %<+% MAMP_elf18_hit_data +
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
                     limits = c(0,100)) +
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Immunogenicity), width = 0.07, offset = 0.12, axis.params = list(line.color = "black"), show.legend = FALSE) +
  scale_fill_manual("Immunogenicity", values = Immunogenicity_colors) 






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

# subset csp22 variant data
MAMP_csp22_hit_data <- subset(filtered_hold_MAMP_seqs, MAMP_Hit == "csp22_consensus")

# add info for ROS outcome 
MAMP_csp22_hit_data["Immunogenicity"] <- "NT"
for (i in 1:nrow(csp22_ros_screen_outcome)){
  MAMP_csp22_hit_data[MAMP_csp22_hit_data$MAMP_Sequence %in% csp22_ros_screen_outcome$Sequence[i],9] <- csp22_ros_screen_outcome$Outcome[i]
} 


# add in new label to match with protein tree tips
MAMP_csp22_hit_data["New_label"] <- paste(MAMP_csp22_hit_data$Protein_Name, MAMP_csp22_hit_data$File_Name, sep = "|")
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[c("New_label", "Percent_Identity", "Genera", "Immunogenicity")]
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[MAMP_csp22_hit_data$New_label %in% CSP_full_protein_tree$tip.label,]


CSP_tree <- ggtree(CSP_full_protein_tree, layout = "circular", ladderize = T, size = 0.2, linetype = 1) %<+% MAMP_csp22_hit_data +
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
  geom_fruit(geom = geom_tile, mapping=aes(fill = Immunogenicity), width = 0.2, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  scale_fill_manual("Immunogenicity", values = Immunogenicity_colors) 



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

