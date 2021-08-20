#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# load all protein sequences with protein annotations of common MAMPs
######################################################################


  All_target_by_annotation <- data.frame("width" = numeric(0), "names" = character(0), "seq" = character(0), 
                                         "Genera" = character(0), "Strain_Name" = character(0), "Filename" = character(0), "Gram" = character(0))
  
  # description to filter gene descriptions for those associated with the MAMPs we are interested in
  key_gene_description <- c("cold-shock","elongation factor Tu","flagellin","cold shock")
  
  for (i in 1:length(load_protein_fasta_files)){
    read_protein_fasta <-  dss2df(Biostrings::readAAStringSet(load_protein_fasta_files[[i]]))
    
    # filter proteins in each fasta for those with the same locus tag as the blast results
    grab_right_protein_seq_blast_results <- read_protein_fasta[grepl(paste(key_gene_description, collapse = "|"), 
                                                                     read_protein_fasta$names),]
    
    # cross reference with metadata to provide genus, strain, filename info
    get_accession_number <- strsplit(load_protein_fasta_files[[i]], "/")[[1]][6]
    get_accession_number <- paste(strsplit(get_accession_number,"_")[[1]][1],
                                  strsplit(get_accession_number,"_")[[1]][2], sep = "_")
    get_strain_info <- subset(datasettable, Assembly_Accession == get_accession_number)
    
    get_accession_info <- cbind(grab_right_protein_seq_blast_results, 
                                rep(get_strain_info$Genera, nrow(grab_right_protein_seq_blast_results)),
                                rep(get_strain_info$Strain_Name, nrow(grab_right_protein_seq_blast_results)),
                                rep(get_strain_info$Filename, nrow(grab_right_protein_seq_blast_results)),
                                rep(get_strain_info$Gram, nrow(grab_right_protein_seq_blast_results)))
    
    All_target_by_annotation <- rbind(All_target_by_annotation, get_accession_info)
  }


  colnames(All_target_by_annotation) <- c("width", "names", "seq", "Genera", "Strain_Name", "Filename", "Gram")
  
  
######################################################################
# Pull out WP tag to filter out proteins which have already been found by blast search
######################################################################


  
  # split All_target_by_annotation to seperate WP tag from rest of protein name
  
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
  correct_blast_df <- data.frame("Sequence" = character(nrow(All_target_by_annotation)), "Percent_Identity" = numeric(nrow(All_target_by_annotation)), "MAMP_Hit" = character(nrow(All_target_by_annotation)))
  pb <- txtProgressBar(min = 0, max = nrow(All_target_by_annotation), style = 3)
  
  
  for (j in 1:nrow(All_target_by_annotation)){
    
    # cold shock protein
    if (grepl(paste(c("cold-shock","cold shock"), collapse = "|"),All_target_by_annotation$names[j]) == TRUE){
      pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'csp22_consensus')
      Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                      All_target_by_annotation$seq[j], type = "global-local", 
                                                                      gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM62)
      correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
      correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
      correct_blast_df$MAMP_Hit[j] <- 'csp22_consensus'
    }
    
    # elongation factor
    if (grepl("elongation factor Tu", All_target_by_annotation$names[j]) == TRUE){
      pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'elf18_consensus')
      Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                      All_target_by_annotation$seq[j], type = "global-local", 
                                                                      gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM62)
      correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
      correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
      correct_blast_df$MAMP_Hit[j] <- 'elf18_consensus'
    }
    
    #flagellin
    if (grepl("flagellin",All_target_by_annotation$names[j]) == TRUE){
      pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flg22_consensus')
      Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                      All_target_by_annotation$seq[j], type = "global-local", 
                                                                      gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM62)
      correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
      correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
      correct_blast_df$MAMP_Hit[j] <- 'flg22_consensus'
    }
    
    setTxtProgressBar(pb, j)
    
  }
  
  close(pb)
  All_target_by_annotation <- cbind(All_target_by_annotation, correct_blast_df)
  
  
  colnames(All_target_by_annotation) <- c("width", "names", "seq", "Genera", "Strain_Name", "Filename", "Gram", "Protein_Name","Sequence","Percent_Identity","MAMP_Hit")
  

