#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


##############################################
# Install packages if needed 
##############################################

#If the required R packages are not found, install them
#If they are installed you can delete these three lines
list_of_packages <- c('rJava','ggnewscale', 'phangorn','tidyverse',
                      'xlsx','ggplot2','ggtree','treeio','tidyr',
                      'ggfortify')
install.packages(list_of_packages)

#Note: If you have issues installing ggtree and treeio from CRAN, 
#      try installing them using BiocManager


#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version="3.9")
BiocManager::install("ggtree")
BiocManager::install("treeio")

##############################################
# load the required packages
##############################################

library(rJava)
library(ggnewscale)
library(phangorn)
library(treeio)
library(ggtree)
library(ggplot2)
library(tidyr)
library(ggfortify)
library(dplyr)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


##############################################
# Load colors - Set tip labels to match other figures
##############################################

#make sure to set path to the same place where the figure 
source("./figure_colors.R")

##datasettable <- datasettable[,c(7,1:6),]
 
if(exists("datasettable") == FALSE){
  source("./tree_tip_data.R")
}

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