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
# Load packages and functions
##############################################

  #make sure to set path to the same place where the figure 
  source("./01_Package_dependencies.R")
  
  source("./02_CommonFunctions.R")
  

##############################################
# Load colors and ggplot theme
##############################################
  
  source("./03_Figure_colors.R")
  
  source("./04_Theme_ggplot.R")
  

##############################################
# Load data and Run through different sections of each script to process data MAMP data
##############################################

  # load data from blast results
  source("./05_Loading_raw_data.R")
  
  
  # we then will processes the BLAST results such that they will be organize in a data table
  source("./06_Process_MAMP_BLAST_results.R")
  
  
  # to be through in our search for microbial MAMPs, we will go back through the annotated genes and pull out the MAMPs that we might have missed in the BLAST search (esspecially for flg22)
  source("./07_Find_MAMPs_by_annotation.R")
  

  # mamps by blast is in "hold_MAMP_seq" and mamps found by annotation are in "All_target_by_annotation"
  source("./08_Combine_blast_and_annotation_results.R")


  # filter for partial proteins and mannual remove a small number of off-targets
  source("./09_Filter_for_partial_proteins.R")


##############################################
# comparing similarity of genomes to filter for clonality
##############################################

  # write out files to run ANI
  source("./10_Parse_genomes_for_ANI_analysis.R")


  #*******************************************
  # go back to methods.md to run ANI analysis
  #*******************************************

  # parsing ANI values and comparing them to MAMP sequences to remove clonal genomes - takes many hours (working on logic to speed up)
  source("./11_ANI_analysis.R")  

  
  # ANI Figure 
  source("./12_ANI_plots.R")


##############################################
# Violin plot Figures - MAMP in comparison to their consensus
##############################################

  # ggplot figures such as violin plots
  source("./13_Comparison_between_MAMP_and_consensus_ggplot_figures.R")




##############################################
# parsing csp22 peptides for classification
##############################################



##############################################
# writing MAMP hits to fasta files for protein tree building
##############################################

  source("./Write_MAMP_hits_to_fasta.R")




##############################################
# make figures
##############################################
  

  
  # Tree's to plot
  source("./Tree_plotting.R")



