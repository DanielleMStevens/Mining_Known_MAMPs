#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# set path to data
######################################################################

#setwd to where repo was cloned and maintained
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
#getSrcDirectory(function(x) {x})



##############################################
# Load colors - Set tip labels to match other figures
##############################################

#make sure to set path to the same place where the figure 
source("./package_dependencies.R")

source("./figure_colors.R")

source("./Theme_ggplot.R")

source("./CommonFunctions.R")

######################################################################
# upload raw data 
######################################################################

# load all files
load_protein_fasta_files <- list.files(path = "./../Protein_fasta/", pattern = ".faa")
load_MAMP_blast_result_files <- list.files(path = './../Protein_fasta/', pattern = ".faa.txt")