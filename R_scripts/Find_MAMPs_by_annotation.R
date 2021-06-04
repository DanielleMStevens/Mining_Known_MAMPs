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

source("./CommonFunctions.R")

######################################################################
# load all protein sequences with protein annotations of common MAMPs
######################################################################


All_target_by_annotation <- data.frame("width" = numeric(0), "names" = character(0), "seq" = character(0))
key_gene_description <- c("cold-shock","elongation factor Tu","flagellin","cold shock")

for (i in 1:length(load_protein_fasta_files)){
  read_protein_fasta <-  dss2df(Biostrings::readAAStringSet(load_protein_fasta_files[[i]]))
  
  # filter proteins in each fasta for those with the same locus tag as the blast results
  grab_right_protein_seq_blast_results <- read_protein_fasta[grepl(paste(key_gene_description, collapse = "|"), 
                                                                   read_protein_fasta$names),]
  
  All_target_by_annotation <- rbind(All_target_by_annotation, grab_right_protein_seq_blast_results)
}


######################################################################
# Pull out WP tag to filter out proteins which have already been found by blast search
######################################################################


# split All_target_by_annotatio  to seperate WP tag from rest of protein name
hold_WP_tag <- list()
for (i in 1:nrow(All_target_by_annotation)){
  hold_WP_tag[[i]] <- strsplit(All_target_by_annotation$names[i], " ", fixed = T)[[1]][1]
  length_of_tag <- length(strsplit(All_target_by_annotation$names[1], " ", fixed = T)[[1]])
  protein_new_name <- paste(strsplit(All_target_by_annotation$names[1], " ", fixed = T)[[1]][2:length_of_tag], collapse = " ")
  #All_target_by_annotation$names[i] <- protein_new_name
}

hold_WP_tag <- as.data.frame(unlist(hold_WP_tag))
All_target_by_annotation <- cbind(All_target_by_annotation, hold_WP_tag)


# pull out mamp sequences to see the likelihood these are the right proteins
correct_blast_df <- data.frame("Sequence" = character(nrow(All_target_by_annotation)), "Percent_Identity" = numeric(nrow(All_target_by_annotation)))

for (j in 1:nrow(All_target_by_annotation)){
  # cold shock protein
  if (grepl(paste(c("cold-shock","cold shock"), collapse = "|"),All_target_by_annotation$names[j]) == TRUE){
    pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'csp22_consensus')
    Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                    All_target_by_annotation$seq[j], type = "global-local", 
                                                                    gapOpening = 100, gapExtension = 100)
    correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
    correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
  }
  # elongation factor
  if (grepl("elongation factor Tu", All_target_by_annotation$names[j]) == TRUE){
    pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'elf18_consensus')
    Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                    All_target_by_annotation$seq[j], type = "global-local", 
                                                                    gapOpening = 100, gapExtension = 100)
    correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
    correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
  }
  #flagellin
  if (grepl("flagellin",All_target_by_annotation$names[j]) == TRUE){
    pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flg22_consensus')
    Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                    All_target_by_annotation$seq[j], type = "global-local", 
                                                                    gapOpening = 100, gapExtension = 100)
    correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
    correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
  }
}


All_target_by_annotation <- cbind(All_target_by_annotation, correct_blast_df)


# filter out hits that have already been found by WP?
filtered_list <- All_target_by_annotation[!All_target_by_annotation$`unlist(hold_WP_tag)`%in% hold_MAMP_seqs$Protein_Name,]




