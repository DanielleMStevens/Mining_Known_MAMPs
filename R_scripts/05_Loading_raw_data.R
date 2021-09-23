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
  load_protein_fasta_files <- list.files(path = "./../../Whole_Genomes_v2/protein_fasta/", pattern = ".faa")
  load_MAMP_blast_result_files <- list.files(path = './../MAMP_blast_outputs/', pattern = ".faa.txt")
  
  # filter
  load_protein_fasta_files <- load_protein_fasta_files[!load_protein_fasta_files %in% load_MAMP_blast_result_files]
  
  
  
  #fix path so can load into process script
  for (i in 1:length(load_protein_fasta_files)){
    load_protein_fasta_files[i] <- paste("./../../Whole_Genomes_v2/protein_fasta/", load_protein_fasta_files[i], sep="")
    load_MAMP_blast_result_files[i] <- paste("./../MAMP_blast_outputs/", load_MAMP_blast_result_files[i], sep="")
  }
  
  
  # load reference file to determine percent identity of mined MAMPs
  load_reference_MAMPs_fasta <- Biostrings::readAAStringSet(filepath = "./../MAMP_database/MAMP_elicitor_list.fasta") #may come with warning messgae..ignore for now
  load_reference_MAMPs_fasta <- dss2df(load_reference_MAMPs_fasta)
  
  
  #########################################################################
  # point to folder with DNA contig files, need to all be formated by Species_Strain.fasta (ex. CC_PF008.fasta)
  #########################################################################
  
  
  #if(exists("datasettable") == FALSE){
    tip_data <- "./../Mining_for_known_MAMPs_genome_accession_info.xlsx" #choose Mining_for_known_MAMPs_genome_accession_info.xlsx
    datasettable <- readxl::read_xlsx(tip_data, col_names = T)
    datasettable <- as.data.frame(datasettable[,1:7])

  #}
