#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


######################################################################
#library packages need to load
######################################################################


#---------Packages for DNA/AA string manipulation---------

#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#BiocManager::install("Biostrings") #if this package won't compile, run the following line on the command line: sudo apt-get install libnlopt-dev
library(Biostrings)
data(BLOSUM62)
data(BLOSUM45)
data(BLOSUM80)
library(pegas)
library(adegenet)

#---------Packages for phylogenic tree loading and manipulation---------

library(ggnewscale)
library(phangorn)
library(treeio)
library(ggtree)
#BiocManager::install("ggtreeExtra")
library(ggtreeExtra)
library(ggbreak)
library(ggdendro)

#---------Packages for general plotting---------

library(ggplot2)
library(extrafont)
library(RColorBrewer)
library(scales)
library(viridis)
library(see)

#---------Packages for text/dataframe loading and manipulation---------

library(readxl)
library(data.table)
library(tidyr)
library(dplyr)
library(reshape2)

#---------Packages for heatmaps and other visuals---------
library(ComplexHeatmap)
library(circlize)
library(ggseqlogo)
#download via devtools
library(devtools)
#library(ggmsa)


#BiocManager::install('XVector')
#library(gghalves)
library(XVector)
library(stringr)
library(readr)
library(plotly)
#library(ggpubr)
#library(rstatix)
library(phylotools)
#library(mat)
#library(rJava)
#library(ggfortify)
library(see)
library(cluster)



##############################################
# Install packages if needed 
##############################################

#If the required R packages are not found, install them
#If they are installed you can delete these three lines
#list_of_packages <- c('rJava','ggnewscale', 'phangorn','tidyverse','xlsx','ggplot2','ggtree','treeio','tidyr','ggfortify')
#install.packages(list_of_packages)

#Note: If you have issues installing ggtree and treeio from CRAN, 
#      try installing them using BiocManager


#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(version="3.9")
#BiocManager::install("ggtree")
#BiocManager::install("treeio")

