#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


#-----------------------Figure 3 Plot------------------------------------------------------------------------------


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



# ---------------------------------------------------------- Figure 3F------------------------------------------------------------------------
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






# ---------------------------------------------------------- Supplemental Figure 6C ------------------------------------------------------------------------

# in order to zoom in on two copies of Streptomyces, we need to A) determine which node to collase (likely the one which seperated immunogenic G+ and G- EFTU)
# sequences. Then replot for supplemental Figure 6C

# label the nodes - this is for visualization only to determine which nodes to expand and collapse
EFTu_tree + geom_text2(aes(subset = !isTip, label = node), hjust = -.3, size = 2) 

# save as 5x5.5 inches - pdf
scaleClade(EFTu_tree, 6451, 0.3) %>% ggtree::collapse(6451, 'min')  %>%  ggtree::collapse(6457, 'min')

# Zoom in on each part of tree - export both as 5 x 3 inches - pdf
ggtree::viewClade(EFTu_tree, 7034) 
ggtree::viewClade(EFTu_tree, 7174) + geom_tiplab(size = 1.5)




# note the code below will add the bootstrap values, however many are too close to the tips.  

#need to make data table to map bootstrap values onto tree
d <- EFTu_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,]

EFTu_tree <- EFTu_tree +
  geom_nodepoint(data=d,aes(label=label), shape = 21, size = 2, fill = "grey35") 







#-----------------------Figure 5 Plot------------------------------------------------------------------------------

##############################################
# load a tree file in newick tree format
##############################################


#CSP_full_protein_tree <- read.tree("./../../Analyses/Protein_alignments_and_trees/cold_shock_protein/full_tree_csp_domain_aligned2.treefile")
CSP_full_protein_tree <- read.tree("./../../Analyses/Protein_alignments_and_trees/cold_shock_protein/full_tree_csp_domain_aligned2.treefile")
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


CSP_tree <- ggtree(CSP_full_protein_tree, layout = "rectangular", ladderize = T, size = 0.2, linetype = 1) %<+% MAMP_csp22_hit_data +
  geom_treescale(x = -10, y = 2, linesize = 1, family = "Arial", offset = 45) +
  geom_tippoint(aes(color = Genera), size = 0.05, alpha = 0.5, show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) +
  

  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.2, offset = 0.1, show.legend = FALSE) +
  scale_fill_manual("Genera", values = Genera_colors) +
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.2, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) +
  
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Immunogenicity), width = 0.2, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  scale_fill_manual("Immunogenicity", values = Immunogenicity_colors) +
  geom_text2(aes(subset = !isTip, label = node), hjust = -.3, size = 0.4)



CSP_tree %>% ggtree::collapse(29123, 'max') %>% 
  #maybe this one
  ggtree::collapse(31126, 'max') %>% 
  #agro clade
  ggtree::collapse(33073, 'max') %>%
  #cspB clade -strep
  ggtree::collapse(32935, 'max') %>% 
  #cspB - clavi
  ggtree::collapse(34260, 'max') %>%
  #cspB -curto
  ggtree::collapse(33971, 'max') %>%
  #cspB -rhodo
  ggtree::collapse(33824, 'max') %>%
  #immo - erwin, pect, dicke
  ggtree::collapse(20136, 'max') %>%
  #agro- immoun
  ggtree::collapse(37798, 'max') %>%
  #pectp - immno
  ggtree::collapse(37356, 'max') %>%
  #agro- long clade
  ggtree::collapse(37796, 'max') 



################

CSP_tree <- ggtree(CSP_full_protein_tree, layout = "rectangular", ladderize = T, size = 0.22, linetype = 1) %<+% MAMP_csp22_hit_data +
  geom_treescale(x = -2, y = 2, linesize = 1, family = "Arial", offset = 45) +
  geom_tippoint(aes(color = Genera), size = 0.03, show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) +
  
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Genera), width = 0.25, offset = 0.1, show.legend = FALSE) +
  scale_fill_manual("Genera", values = Genera_colors) +
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.25, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  scale_fill_viridis(option="magma", name = "Percent AA\nSimilarity", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) +
  
  
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Immunogenicity), width = 0.25, offset = 0.1, axis.params = list(line.color = "black"), show.legend = FALSE) +
  scale_fill_manual("Immunogenicity", values = Immunogenicity_colors) 
  #geom_text2(aes(subset = !isTip, label = node), hjust = -.3, size = 0.4)


### 4x8.8


# ---------------------------------------------------------- Supplemental Figure 6C ------------------------------------------------------------------------



MPMI_tree <- ggtree::ggtree(core_gene_phylo, ladderize = T, size = 0.4, linetype = 1) %<+% map_data +
  geom_tippoint(aes(color = Genera), size = 0.3, show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) +
  layout_fan(angle=90) +
  geom_text2(aes(subset = !isTip, label = node), hjust = 1, size = 0.4)


MPMI_tree <- ggtree::ggtree(core_gene_phylo, ladderize = T, size = 0.4, linetype = 1, layout = 'circular') %<+% map_data +
  geom_tippoint(aes(color = Genera), size = 0.3, show.legend = FALSE) +
  scale_color_manual("Genera", values = Genera_colors) 
  #geom_text2(aes(subset = !isTip, label = node), hjust = 1, size = 0.4)

ggtree::viewClade(MPMI_tree, 6813) 



# need to make data table to map bootstrap values onto tree
d <- csp22_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 70,]


