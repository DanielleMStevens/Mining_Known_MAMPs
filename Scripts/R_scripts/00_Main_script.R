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
  library(rstudioapi)
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
  
  
  # to be through in our search for microbial MAMPs, we will go back through the annotated genes and pull out the MAMPs 
  # that we might have missed in the BLAST search (esspecially for flg22)
  source("./07_Find_MAMPs_by_annotation.R")
  

  # mamps by blast is in "hold_MAMP_seq" and mamps found by annotation are in "All_target_by_annotation"
  source("./08_Combine_blast_and_annotation_results.R")


  # filter for partial proteins and mannual remove a small number of off-targets
  source("./09_Filter_for_partial_proteins.R")



##############################################
# comparing similarity of genomes to filter for clonality
##############################################

  # write out files to run ANI - run ONCE to generate file to run through fastANI
  source("./10_Parse_genomes_for_ANI_analysis.R")


  #*******************************************
  # go back to methods.md to run ANI analysis
  #*******************************************

  # parsing ANI values and comparing them to MAMP sequences to remove clonal genomes
  source("./11_ANI_analysis.R")  
  
  
  # finalizing MAMP list
  source("./12_Finalize_MAMP_list_post_ANI.R")
  
  
  # ANI Figure  - takes a long time to run so only run once.
  # Part of Supplemental Figure 1
  source("./13_ANI_plots.R")
  

##############################################
# Dyanamic changes of MAMPs from diverse bacteria - Figure 1 & 2
##############################################
  
  #*******************************************
  # go back to methods.md to run phylogenetic analysis
  #*******************************************
  
  # phylogenetic tree of diverse bacteria with MAMP eptitope number plotted on
  source("./14_Core_tree_with_MAMP_number_plotted.R")

  
  # ggplot figures such as violin plots
  # includes Figure 1C, D, E; as well as part of Supplemental Figure 1
  source("./15_Comparison_between_MAMP_and_consensus_ggplot_figures.R")
  
  
  # ggplot figure for nlp20 - Supplemental Figure 1
  source("./16_nlp20_plots.R")
  
  
  # A deeper dive on eptiope variation at the genera level
  source("./17_Compare_mamps_to_consensus_by_genera_figures.R")
  
  
  # Assessing variation within each MAMP epitope and in a positional manner 
  ### this takes a long time to complete -1-1.5 hours
  source("./18_Assessing_MAMP_variation.R")
  
  
  # how many MAMPs are there, the variation that extists, and how many/offen plot
  # includes Supplemental Figure 4
  source("./19_Abundance_of_MAMPs.R")


##############################################
# Diverse epitopes and their impact on plant immmune perception - Figure 2 and 3
##############################################

  # Creating phylogenic tree with immunogenicity data (via ROS assays)
  source("./20_ROS_Screen_tree.R")
  
##############################################
# writing MAMP hits to fasta files for protein tree building
##############################################
  
  # Re-pull whole protein sequences for MAMP hits and write to fasta file
  source("./21_Write_MAMP_hits_to_fasta.R")
  
  #*******************************************
  # go back to methods.md to run phylogenetic analysis - build protein trees
  #*******************************************
  
  # plot whole protein trees 
  source("./22_Plot_protein_tree")
  
  
##############################################
# parsing csp22 peptides for classification
##############################################

  # to classify CSPs, first we will assess those within each genera
  # Pull CSPs from each genera and write to a fasta file
  source("./23_by_genera_csps_to_fasta.R")

  # combining mmseq2, meme, and protein tree analyses to determine number of csp types
  source("./24_CSP_classification.R")

  # comparing csp types across genera
  source("./25_ortho_csp_cluster.R")




