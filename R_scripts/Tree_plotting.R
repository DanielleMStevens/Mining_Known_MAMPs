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
library(tidyverse)
library(tidyr)
library(ggfortify)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


##############################################
# Load colors - Set tip labels to match other figures
##############################################

#make sure to set path to the same place where the figure 
source("./figure_colors.R")

datasettable <- datasettable[,c(7,1:6),]

##############################################
# load a tree file in newick tree format
##############################################


csp22_full_protein_tree <- read.tree("./../Protein_alignments_and_trees/csp22_full_length_alignment.treefile")
csp22_full_protein_tree <- phangorn::midpoint(csp22_full_protein_tree, node.labels='label')



for (i in 1:length(csp22_full_protein_tree$tip.label)){
  csp22_full_protein_tree$tip.label[i] <- strsplit(csp22_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5]
}

test <- datasettable[datasettable$Filename %in% csp22_full_protein_tree$tip.label,]

csp22_tree <- ggtree(csp22_full_protein_tree, layout = 'circular', ladderize = T, size = 0.3, linetype = 1) %<+% test +
  geom_tippoint(aes(color = Genera), size = 0.8) +
  scale_color_manual("Genera", values = Genera_colors) +
  theme(legend.position = "none")


# need to make data table to map bootstrap values onto tree
d <- csp22_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 70,]

csp22_tree <- csp22_tree + 
  geom_nodepoint(data = d, aes(label = label),size = 0.6, alpha = 0.6)



#######################################################################
# fix names on tip labels
#######################################################################



elf18_full_protein_tree <- read.tree("./../Protein_alignments_and_trees/EfTu/elf18_full_length_alignment.treefile")
elf18_full_protein_tree <- phangorn::midpoint(elf18_full_protein_tree, node.labels='label')



for (i in 1:length(elf18_full_protein_tree$tip.label)){
  elf18_full_protein_tree$tip.label[i] <- strsplit(elf18_full_protein_tree$tip.label[i], "|", fixed = T)[[1]][5]
}

test <- datasettable[datasettable$Filename %in% elf18_full_protein_tree$tip.label,]

elf18_tree <- ggtree(elf18_full_protein_tree, layout = 'circular', ladderize = T, size = 0.3, linetype = 1) %<+% test +
  geom_tippoint(aes(color = Genera), size = 0.8) +
  scale_color_manual("Genera", values = Genera_colors) +
  theme(legend.position = "none")


# need to make data table to map bootstrap values onto tree
d <- csp22_tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 70,]

csp22_tree <- csp22_tree + 
  geom_nodepoint(data = d, aes(label = label),size = 0.6, alpha = 0.6)

