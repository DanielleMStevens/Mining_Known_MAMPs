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



##############################################
# Load colors - Set tip labels to match other figures
##############################################

  #make sure to set path to the same place where the figure 
  source("./package_dependencies.R")
  

  source("./CommonFunctions.R")
  


##############################################
# Run through different sections of each script to process data
##############################################

  # load data from blast results
  source("./Loading_raw_data.R")
  
  
  # we then will processes the BLAST results such that they will be organize in a data table
  source("./Process_MAMP_BLAST_results.R")
  
  
  # to be through in our search for microbial MAMPs, we will go back through the annotated genes and pull out
  # the MAMPs that we might have missed in the BLAST search (esspecially for flg22)
  source("./Find_MAMPs_by_annotation.R")
  

  # mamps by blast is in "hold_MAMP_seq" and mamps found by annotation are in "All_target_by_annotation"
  source("./combine_blast_and_annotation_results.R")

  
##############################################
# parsing csp22 peptides for classification
##############################################



##############################################
# writing MAMP hits to fasta files for protein tree building
##############################################

  source("./Write_MAMP_hits_to_fasta.R")


##############################################
# Load colors - Set tip labels to match other figures
##############################################

  source("./figure_colors.R")
  
  source("./Theme_ggplot.R")

##############################################
# make figures
##############################################
  
  # ANI FIgure - will be used for tree ploting figure
  source("./ANI_analysis.R")
  
  # Tree's to plot
  source("./Tree_plotting.R")

  # ggplot figures such as violin plots
  source("./ggplot_figures.R")




