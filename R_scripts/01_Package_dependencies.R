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

#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#BiocManager::install("Biostrings") #if this package won't compile, run the following line on the command line: sudo apt-get install libnlopt-dev
#BiocManager::install('XVector')
#BiocManager::install("ggtreeExtra")


library(ggplot2)
library(extrafont)
library(devtools)
library(readxl)
library(reshape2)
library(rstudioapi)
library(Biostrings)
#library(ggmsa)
#library(ggbeeswarm)
#library(gghalves)
library(XVector)
library(RColorBrewer)
library(devtools)
library(readxl)
library(stringr)
library(readr)
library(rstudioapi)
library(plotly)
#library(ggpubr)
#library(rstatix)
library(phylotools)
library(scales)
#library(mat)
#library(rJava)
library(ggnewscale)
library(phangorn)
library(treeio)
library(ggtree)
library(ggplot2)
library(tidyr)
#library(ggfortify)
library(dplyr)
library(ggtreeExtra)
library(viridis)
library(see)


data(BLOSUM62)
data("BLOSUM45")
data("BLOSUM80")

library(ComplexHeatmap)
library(reshape2)
library(stringr)
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

##############################################
# load the required packages
##############################################




#library(Biostrings)