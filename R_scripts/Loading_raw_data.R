#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# upload raw data 
######################################################################

  # load all files
  load_protein_fasta_files <- list.files(path = "./../../Whole_Genomes/protein_fasta/", pattern = ".faa")
  load_MAMP_blast_result_files <- list.files(path = './../MAMP_blast_outputs/', pattern = ".faa.txt")
  
  # filter
  load_protein_fasta_files <- load_protein_fasta_files[!load_protein_fasta_files %in% load_MAMP_blast_result_files]
  
  
  
  #fix path so can load into process script
  for (i in 1:length(load_protein_fasta_files)){
    load_protein_fasta_files[i] <- paste("./../../Whole_Genomes/protein_fasta/", load_protein_fasta_files[i], sep="")
    load_MAMP_blast_result_files[i] <- paste("./../MAMP_blast_outputs/", load_MAMP_blast_result_files[i], sep="")
  }
  
  
  # load reference file to determine percent identity of mined MAMPs
  load_reference_MAMPs_fasta <- Biostrings::readAAStringSet(filepath = "./../MAMP_database/MAMP_elicitor_list.fasta") #may come with warning messgae..ignore for now
  load_reference_MAMPs_fasta <- dss2df(load_reference_MAMPs_fasta)
  
