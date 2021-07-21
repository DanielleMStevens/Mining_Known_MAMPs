#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------




##############################################
# Load colors - Set tip labels to match other figures
##############################################

##datasettable <- datasettable[,c(7,1:6),]
 
#if(exists("datasettable") == FALSE){
# source("./tree_tip_data.R")
#}


##############################################
# ANI based phylogeny
##############################################

row_dend = hclust(dist(empty_ANI_matrix))

# change out label names
for (i in 1:length(row_dend$labels)){
  accession_label <- paste(strsplit(row_dend$labels[[i]], "_")[[1]][1],
                           strsplit(row_dend$labels[[i]], "_")[[1]][2],
                           sep = "_")
  row_dend$labels[[i]] <- datasettable[datasettable$Assembly_Accession %in% accession_label,7]
}

map_data <- datasettable[,c(7,1:6)]

csp22_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "csp22_consensus")
map_data <- cbind(map_data, "csp22_consensus" = csp22_copy_number[match(map_data[,1], csp22_copy_number[,2]),3])

elf18_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "elf18_consensus")
map_data <- cbind(map_data, "elf18_consensus" = elf18_copy_number[match(map_data[,1], elf18_copy_number[,2]),3])

flg22_copy_number <- subset(hold_copy_number, hold_copy_number$MAMP_Hit == "flg22_consensus")
map_data <- cbind(map_data, "flg22_consensus" = flg22_copy_number[match(map_data[,1], flg22_copy_number[,2]),3])
map_data[is.na(map_data$flg22_consensus),10] <- 0



ggtree::ggtree(as.dendrogram(row_dend), layout = 'circular', ladderize = T, 
               size = 0.35, linetype = 1) %<+% map_data +
  geom_fruit(geom = geom_tile, mapping = aes(color = Genera, fill = Genera), width = 30,
             offset = 0.02, axis.params = list(line.color = "black")) +
  scale_color_manual("Genera", values = Genera_colors) +
  scale_fill_manual("Genera", values = Genera_colors) +
  ggnewscale::new_scale_fill() +

  geom_fruit(geom = geom_tile, mapping = aes(fill = csp22_consensus), width = 35,
             offset = 0.05) +
  geom_fruit(geom = geom_tile, mapping = aes(fill = elf18_consensus), width = 35,
             offset = 0.05) + 
  geom_fruit(geom = geom_tile, mapping = aes(fill = flg22_consensus), width = 35,
             offset = 0.05) +
  scale_fill_gradient(low = "white", high = "black", breaks = c(0,2,4,6,8,10,12,14), guide = "legend")





##############################################
# core gene phylogeny
##############################################


core_gene_phylo <- read.tree(file.choose())
core_gene_phylo <- phangorn::midpoint(core_gene_phylo, node.labels='label')

for (i in 1:length(core_gene_phylo$tip.label)){
  accession_label <- paste(strsplit(core_gene_phylo$tip.label[[i]], "_")[[1]][1],
                           strsplit(core_gene_phylo$tip.label[[i]], "_")[[1]][2],
                           sep = "_")
  core_gene_phylo$tip.label[[i]] <- datasettable[datasettable$Assembly_Accession %in% accession_label,1]
}
  
#layout="circular", ladderize = T, size = 0.35, linetype = 1
ggtree(core_gene_phylo)

##############################################
# load a tree file in newick tree format
##############################################


csp22_full_protein_tree <- read.tree("./../Protein_alignments_and_trees/cold_shock_protein/csp22_full_length_alignment.treefile")
csp22_full_protein_tree <- phangorn::midpoint(csp22_full_protein_tree, node.labels='label')


for (i in 1:length(csp22_full_protein_tree$tip.label)){
  csp22_full_protein_tree$tip.label[i] <- paste(strsplit(csp22_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                                strsplit(csp22_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                                sep = "|")
}


MAMP_csp22_hit_data <- subset(hold_MAMP_seqs, MAMP_Hit == "csp22_consensus")
MAMP_csp22_hit_data["New_label"] <- paste(MAMP_csp22_hit_data$Protein_Name, MAMP_csp22_hit_data$File_Name, sep = "|")
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[,c(11,3,7)]
MAMP_csp22_hit_data <- MAMP_csp22_hit_data[MAMP_csp22_hit_data$New_label %in% csp22_full_protein_tree$tip.label,]


csp22_tree <- ggtree(csp22_full_protein_tree, layout="circular", ladderize = T, size = 0.35, linetype = 1) %<+% MAMP_csp22_hit_data +
  geom_tippoint(aes(color = Genera), size = 0.5) +
  geom_tiplab(size = 0.3, align = T, linesize = 0.015)+
  scale_color_manual("Genera", values = Genera_colors) +
  theme(legend.position = "none") +
  geom_treescale(fontsize = 2, linesize = 0.3, x = 1, y = 0) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.12,
             offset = 0.08, axis.params = list(line.color = "black")) +
  scale_fill_viridis(option="magma", name = "Percent AA Similarity\n", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) 

csp22_tree <- csp22_tree + ggplot2::theme(legend.position = c(0.02, 0.3), legend.title = element_text("Genera", family = "Arial" ,
                                                                               face = "bold", color = "black", size = 12),
                   legend.text = element_text(family = "Arial", color = "black", size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))




csp22_tree


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
# plotting protein tree for EF-Tu
#######################################################################



elf18_full_protein_tree <- read.tree("./../Protein_alignments_and_trees/EfTu/elf18_full_length_alignment.treefile")
elf18_full_protein_tree <- phangorn::midpoint(elf18_full_protein_tree, node.labels='label')



for (i in 1:length(elf18_full_protein_tree$tip.label)){
  elf18_full_protein_tree$tip.label[i] <- paste(strsplit(elf18_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                                strsplit(elf18_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                                sep = "|")
}


MAMP_elf18_hit_data <- subset(hold_MAMP_seqs, MAMP_Hit == "elf18_consensus")
MAMP_elf18_hit_data["New_label"] <- paste(MAMP_elf18_hit_data$Protein_Name, MAMP_elf18_hit_data$File_Name, sep = "|")
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[,c(11,3,7)]
MAMP_elf18_hit_data <- MAMP_elf18_hit_data[MAMP_elf18_hit_data$New_label %in% elf18_full_protein_tree$tip.label,]


elf18_tree <- ggtree(elf18_full_protein_tree, layout="circular", ladderize = T, size = 0.35, linetype = 1) %<+% MAMP_elf18_hit_data +
  geom_tippoint(aes(color = Genera), size = 0.25) +
  scale_color_manual("Genera", values = Genera_colors) +
  geom_tiplab(size = 0.3, align = T, linesize = 0.1)+
  theme(legend.position = "none") +
  geom_treescale(fontsize = 2, linesize = 0.3, x = 2, y = 0) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.12,
             offset = 0.08, axis.params = list(line.color = "black")) +
  scale_fill_viridis(option="magma", name = "Percent AA Similarity\n", 
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c(0, 20, 40, 60, 80, 100),
                     limits = c(0,100)) 
  
elf18_tree + theme(legend.position = c(0.05, 0.3), legend.title = element_text("Genera", family = "Arial" ,
                                                                              face = "bold", color = "black", size = 12),
                   legend.text = element_text(family = "Arial", color = "black", size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))




# need to make data table to map bootstrap values onto tree
d <- elf18_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 70,]


#######################################################################
# plotting protein tree for Flagellin
#######################################################################


flg22_full_protein_tree <- read.tree("./../Protein_alignments_and_trees/Flagellin/flg22_full_length_alignment.treefile")
flg22_full_protein_tree <- phangorn::midpoint(flg22_full_protein_tree, node.labels='label')



for (i in 1:length(flg22_full_protein_tree$tip.label)){
  flg22_full_protein_tree$tip.label[i] <- paste(strsplit(flg22_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                                strsplit(flg22_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                                sep = "|")
}


MAMP_flg22_hit_data <- subset(hold_MAMP_seqs, MAMP_Hit == "flg22_consensus")
MAMP_flg22_hit_data["New_label"] <- paste(MAMP_flg22_hit_data$Protein_Name, MAMP_flg22_hit_data$File_Name, sep = "|")
MAMP_flg22_hit_data <- MAMP_flg22_hit_data[,c(11,3,7)]
MAMP_flg22_hit_data <- MAMP_flg22_hit_data[MAMP_flg22_hit_data$New_label %in% flg22_full_protein_tree$tip.label,]


flg22_tree <- ggtree(flg22_full_protein_tree, layout="circular", ladderize = T, size = 0.35, linetype = 1) %<+% MAMP_flg22_hit_data +
  geom_tippoint(aes(color = Genera), size = 0.5) +
  scale_color_manual("Genera", values = Genera_colors) +
  theme(legend.position = "none") +
  geom_treescale(fontsize = 2, linesize = 0.3, x = 1.5, y = 0) +
  ggnewscale::new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping=aes(fill = Percent_Identity), width = 0.12,
             offset = 0.08, axis.params = list(line.color = "black")) +
  scale_fill_viridis(option="magma", name = "Percent AA Similarity\n", 
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





############# for rectangular treee
#max_value_edge <- max(elf18_full_protein_tree$edge.length)*2.5
theme_tree2() + #this next three lines are a hacky way to "rescale" the tree 
  xlim(0, max_value_edge)+
  theme_tree() +