#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/08/2023
# Script Purpose: Finds all MAMPs by local alignment to annotated genes
# Inputs: Protein Fasta files
# Outputs: N/A, temporary data table which holds all the MAMP sequence (plus additonal info)
#-----------------------------------------------------------------------------------------------



######################################################################
# load all protein sequences with protein annotations of common MAMPs
######################################################################

  # empty table to fill in regarding the MAMP, the Gene, sequence and information on the genome their dervied from
  All_target_by_annotation <- data.frame("width" = numeric(0), "names" = character(0), "seq" = character(0), 
                                         "Genera" = character(0), "Strain_Name" = character(0), "Filename" = character(0), "Gram" = character(0))
  
  # description to filter gene descriptions for those associated with the MAMPs we are interested in
  key_gene_description <- c("cold-shock","elongation factor Tu","flagellin","cold shock","NPP1","RNA chaperone/antiterminator","transcription antiterminator/RNA stability")
  
  for (i in 1:length(load_protein_fasta_files)){
    read_protein_fasta <-  dss2df(Biostrings::readAAStringSet(load_protein_fasta_files[[i]]))
    
    # filter proteins in each fasta for those with the same locus tag as the blast results
    grab_right_protein_seq_blast_results <- read_protein_fasta[grepl(paste(key_gene_description, collapse = "|"), 
                                                                     read_protein_fasta$names),]
    
    # cross reference with metadata to provide genus, strain, filename info
    get_accession_number <- Biostrings::strsplit(load_protein_fasta_files[[i]], "/")[[1]][7]
    get_accession_number <- paste(strsplit(get_accession_number,"_")[[1]][1],
                                  strsplit(get_accession_number,"_")[[1]][2], sep = "_")
    get_strain_info <- subset(datasettable, Assembly_Accession == get_accession_number)
    
    get_accession_info <- cbind(grab_right_protein_seq_blast_results, 
                                rep(get_strain_info$Genera, nrow(grab_right_protein_seq_blast_results)),
                                rep(get_strain_info$Strain_Name, nrow(grab_right_protein_seq_blast_results)),
                                rep(get_strain_info$Filename, nrow(grab_right_protein_seq_blast_results)),
                                rep(get_strain_info$Gram, nrow(grab_right_protein_seq_blast_results)))
    
    # add all data into empty table declared above
    All_target_by_annotation <- rbind(All_target_by_annotation, get_accession_info)
  }


  colnames(All_target_by_annotation) <- c("width", "names", "seq", "Genera", "Strain_Name", "Filename", "Gram")
  
  
######################################################################
# Pull out WP tag to filter out proteins which have already been found by blast search
# This allows us to pull MAMP seqs not found by blastp
######################################################################
  
  # split All_target_by_annotation to seperate WP tag from rest of protein name
  hold_WP_tag <- list()
  
  for (i in 1:nrow(All_target_by_annotation)){
    hold_WP_tag[[i]] <- Biostrings::strsplit(All_target_by_annotation$names[i], " ", fixed = T)[[1]][1]
    length_of_tag <- length(strsplit(All_target_by_annotation$names[1], " ", fixed = T)[[1]])
    protein_new_name <- paste(strsplit(All_target_by_annotation$names[1], " ", fixed = T)[[1]][2:length_of_tag], collapse = " ")
    # All_target_by_annotation$names[i] <- protein_new_name
  }
  
  hold_WP_tag <- as.data.frame(unlist(hold_WP_tag))
  All_target_by_annotation <- cbind(All_target_by_annotation, hold_WP_tag)
  

    
  # pull out mamp sequences to see the likelihood these are the right proteins
  correct_blast_df <- data.frame("Sequence" = character(nrow(All_target_by_annotation)), 
                                 "Percent_Identity" = numeric(nrow(All_target_by_annotation)), 
                                 "MAMP_Hit" = character(nrow(All_target_by_annotation)))
  
  correct_blast_df_flgII_28 <- data.frame("width" = numeric(0), "names" = character(0), "seq" = character(0), "Genera" = character(0),
                                          "Strain_Name" = character(0), "Filename" = character(0), "Gram" = character(0), "Protein_Name" = character(0),
                                          "Sequence" = character(0), "Percent_Identity" = numeric(0), "MAMP_Hit" = character(0))
  
  
  
  for (j in 1:nrow(All_target_by_annotation)){
    
    # -----------------------cold shock protein-----------------------
    if (grepl(paste(c("cold-shock","cold shock","RNA chaperone/antiterminator","transcription antiterminator/RNA stability"), collapse = "|"), All_target_by_annotation$names[j]) == TRUE){
      
      # if csp protein is less than 50 AA long, it's very likely a partial CDS and will be removed
      if (All_target_by_annotation$width[j] < 50){
        next
      }
      
      # if csp protein is longer than 50 AA long, its more likely to be a full length sequence
      if (All_target_by_annotation$width[j] > 50){
        pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'csp22_consensus')
        Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment("KGFGF", 
                                                                        All_target_by_annotation$seq[j], type = "global-local", 
                                                                        gapOpening = 300, substitutionMatrix = BLOSUM80, scoreOnly = T)
        
        # ignore if score is low (i.e. doesn't carry crtical motif)
        if (Alignment_between_MAMP_and_Ref < 20){
          next
        }
        
        # if score high enough, find mamp via local alignment to consensus
        if (Alignment_between_MAMP_and_Ref > 20){
        
          Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                          All_target_by_annotation$seq[j], type = "global-local", 
                                                                          gapOpening = 300, substitutionMatrix = BLOSUM45)
          # ignore - for degugging purposes
          # print(Alignment_between_MAMP_and_Ref)
          # print(paste(All_target_by_annotation$`unlist(hold_WP_tag)`[j], All_target_by_annotation$seq[j], sep = "|"))
          correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
          correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
          correct_blast_df$MAMP_Hit[j] <- 'csp22_consensus'
        }
      }
    }
    
    # -----------------------elongation factor-----------------------
    if (grepl("elongation factor Tu", All_target_by_annotation$names[j]) == TRUE){
      
      if (grepl("partial", All_target_by_annotation$names[j]) == TRUE){
        next
      }
      if (grepl("partial", All_target_by_annotation$names[j]) == FALSE){
        
        pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'elf18_consensus')
        Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                        All_target_by_annotation$seq[j], type = "global-local", 
                                                                        gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM45)
        if (Alignment_between_MAMP_and_Ref@score < 50){
          next
        }
        if (Alignment_between_MAMP_and_Ref@score > 50){
          
          # ignore - for degugging purposes
          # print(Alignment_between_MAMP_and_Ref)
          # print(All_target_by_annotation$`unlist(hold_WP_tag)`[j])
          correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
          correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
          correct_blast_df$MAMP_Hit[j] <- 'elf18_consensus'
        }
      }
    }
    
    # -----------------------flagellin-----------------------
    if (grepl("flagellin",All_target_by_annotation$names[j]) == TRUE){
      
      if (grepl("partial", All_target_by_annotation$names[j]) == TRUE){
        next
      }
      if (grepl("partial", All_target_by_annotation$names[j]) == FALSE){
      
        pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flg22_consensus')
        Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                        All_target_by_annotation$seq[j], type = "global-local", 
                                                                        gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM45)
        
        if (Alignment_between_MAMP_and_Ref@score > 40){

        
          # print(Alignment_between_MAMP_and_Ref)
          # print(All_target_by_annotation$`unlist(hold_WP_tag)`[j])
          correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
          correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
          correct_blast_df$MAMP_Hit[j] <- 'flg22_consensus'
        }
        
        pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flgII-28')
        Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                        All_target_by_annotation$seq[j], type = "global-local", 
                                                                        gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM62)
        if (Alignment_between_MAMP_and_Ref@score > 40){
          # print(Alignment_between_MAMP_and_Ref)
          # print(All_target_by_annotation$`unlist(hold_WP_tag)`[j])
          
          if (identical(as.character(Alignment_between_MAMP_and_Ref@subject), character(0)) == FALSE){
            temp_df <- data.frame(All_target_by_annotation[j,],
                                  "Sequence" =  as.character(Alignment_between_MAMP_and_Ref@subject),
                                  "Percent_Identity" = Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1"),
                                  "MAMP_Hit" = 'flgII.28')
            colnames(temp_df) <- colnames(correct_blast_df_flgII_28)
            correct_blast_df_flgII_28 <- rbind(correct_blast_df_flgII_28, temp_df)
          }
        }
      }
    }
    
    # -----------------------NPP1 -  necrosis-----------------------
    if (grepl("NPP1", All_target_by_annotation$names[j]) == TRUE){
      pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'nlp20_consensus')
      Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                      All_target_by_annotation$seq[j], type = "global-local", 
                                                                      gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM62)
      
      # print(Alignment_between_MAMP_and_Ref)
      # print(All_target_by_annotation$`unlist(hold_WP_tag)`[j])
      correct_blast_df$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
      correct_blast_df$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
      correct_blast_df$MAMP_Hit[j] <- 'nlp20_consensus'
    }
    
    setTxtProgressBar(pb, j)
    
  }
  
  
######################################################################
# merge and clean-up
######################################################################
  
  
  All_target_by_annotation <- cbind(All_target_by_annotation, correct_blast_df)
  colnames(All_target_by_annotation) <- c("width", "names", "seq", "Genera", "Strain_Name", "Filename", "Gram", "Protein_Name","Sequence","Percent_Identity","MAMP_Hit")
  
  All_target_by_annotation <- rbind(All_target_by_annotation, correct_blast_df_flgII_28)
  

  rm(correct_blast_df)
  rm(correct_blast_df_flgII_28)
  
  #remove hits which were 0'd since they are so low in value, they likely aren't real
  filter_list <- list()
  for (i in 1:nrow(All_target_by_annotation)){
    if (All_target_by_annotation$Percent_Identity[i] == 0.00000){
      filter_list[[i]] <- FALSE
    }
    if (All_target_by_annotation$Percent_Identity[i] != 0.00000){
      filter_list[[i]] <- TRUE
    }
  }
  
  All_target_by_annotation <- All_target_by_annotation[unlist(filter_list),]
  rm(filter_list)
  
  
  
  
  

