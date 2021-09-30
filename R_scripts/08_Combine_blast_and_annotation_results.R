#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# merge two data bases
######################################################################

# parse any annotaiton hits that have already been find by blast
annotation_list <- paste(All_target_by_annotation$Genera, All_target_by_annotation$Strain_Name,
                         All_target_by_annotation$Protein_Name, All_target_by_annotation$MAMP_Hit, sep = "|")

blast_list <- paste(hold_MAMP_seqs$Genera, hold_MAMP_seqs$Strain_Name, hold_MAMP_seqs$Protein_Name, hold_MAMP_seqs$MAMP_Hit, sep = "|")

filter_list <- All_target_by_annotation[!annotation_list %in% blast_list,]



# align both databases so they can be merged
hold_MAMP_seqs <- hold_MAMP_seqs[,c(1,2,3,6,7,8,9,10)]
filter_list <- filter_list[c(8,11,10,9,4,5,6,7)]
rownames(filter_list) <- NULL
colnames(filter_list) <- colnames(hold_MAMP_seqs)
hold_MAMP_seqs <- rbind(hold_MAMP_seqs, filter_list)  

rm(filter_list)

#hold_MAMP_seqs <- subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit != "nlp20_consensus")


######################################################################
# correct for off targets
######################################################################

    # can reuse annotation database to use full-length protein by WP tag to adjust for hit position
    
    ######################################################################
    # load all protein sequences with protein annotations of common MAMPs
    ######################################################################
    
    
    All_target_by_annotation <- data.frame("width" = numeric(0), "names" = character(0), "seq" = character(0), 
                                           "Genera" = character(0), "Strain_Name" = character(0), "Filename" = character(0), "Gram" = character(0))
    
    # description to filter gene descriptions for those associated with the MAMPs we are interested in
    key_gene_description <- c("cold-shock","elongation factor Tu","flagellin","cold shock","NPP1")
    
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
    
    
    # fix flgII-28 name scheme in hold_MAMP seq
    for (i in 1:nrow(hold_MAMP_seqs)){
      if (hold_MAMP_seqs$MAMP_Hit[i] == "flgII.28"){
        hold_MAMP_seqs$MAMP_Hit[i] <- "flgII-28"
      }
    }
    
    
    
    # readjust hit sequence if value seems off
    filter_list_2 <- list()
    for (i in 1:nrow(hold_MAMP_seqs)){
      
      if (hold_MAMP_seqs$Percent_Identity[i] > 40){
        filter_list_2[[i]] <- TRUE
      }
      
      if (hold_MAMP_seqs$Percent_Identity[i] <= 40){
        
        # cold shock protein
        if (hold_MAMP_seqs$MAMP_Hit[i] == "csp22_consensus"){
          
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'csp22_consensus')
          hold_MAMP_seqs$MAMP_Sequence[i] <- str_replace_all(hold_MAMP_seqs$MAMP_Sequence[i], "-", "")
          Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment("KGFGF", 
                                                                          hold_MAMP_seqs$MAMP_Sequence[i], type = "global-local", 
                                                                          gapOpening = 300, substitutionMatrix = BLOSUM80, scoreOnly = TRUE)
          

          
          if (Alignment_between_MAMP_and_Ref < 20){
            filter_list_2[[i]] <- FALSE
          }
          
          if (Alignment_between_MAMP_and_Ref > 20){
            
            filter1 <- All_target_by_annotation[All_target_by_annotation$Genera %in% hold_MAMP_seqs$Genera[i],]
            filter2 <- filter1[filter1$Strain_Name %in% hold_MAMP_seqs$Strain_Name[i],]
            filter3 <- filter2[filter2$`unlist(hold_WP_tag)` %in% hold_MAMP_seqs$Protein_Name[i],]

            
            if (nrow(filter3) == 0){
              filter_list_2[[i]] <- FALSE
            }
            
            if (nrow(filter3) != 0){
              Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                              filter3$seq[1], type = "global-local", 
                                                                              gapOpening = 300, substitutionMatrix = BLOSUM45)

              percent_identity_realignment <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")

              
              if (hold_MAMP_seqs$Percent_Identity[i] >= percent_identity_realignment){
                filter_list_2[[i]] <- TRUE
              }
              if (hold_MAMP_seqs$Percent_Identity[i] < percent_identity_realignment){
                hold_MAMP_seqs$MAMP_Sequence[i] <- as.character(Alignment_between_MAMP_and_Ref@subject)
                hold_MAMP_seqs$Percent_Identity[i] <- percent_identity_realignment
                filter_list_2[[i]] <- TRUE
              }
            }
          }
        }
        
        
        # elongation factor
        if (hold_MAMP_seqs$MAMP_Hit[i] == "elf18_consensus"){
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'elf18_consensus')
          hold_MAMP_seqs$MAMP_Sequence[i] <- str_replace_all(hold_MAMP_seqs$MAMP_Sequence[i], "-", "")
          Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment("NXGTXG", 
                                                                          hold_MAMP_seqs$MAMP_Sequence[i], type = "global-local", 
                                                                          gapOpening = 300, substitutionMatrix = BLOSUM80, scoreOnly = TRUE)
          #print(Alignment_between_MAMP_and_Ref)
          #print(hold_MAMP_seqs[i,])
          
          if (Alignment_between_MAMP_and_Ref < 20){
            filter_list_2[[i]] <- FALSE
          }
          
          if (Alignment_between_MAMP_and_Ref > 20){
            
            filter1 <- All_target_by_annotation[All_target_by_annotation$Genera %in% hold_MAMP_seqs$Genera[i],]
            filter2 <- filter1[filter1$Strain_Name %in% hold_MAMP_seqs$Strain_Name[i],]
            filter3 <- filter2[filter2$`unlist(hold_WP_tag)` %in% hold_MAMP_seqs$Protein_Name[i],]
            
            
            if (nrow(filter3) == 0){
              filter_list_2[[i]] <- FALSE
            }
            
            if (nrow(filter3) != 0){
              Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                              filter3$seq[1], type = "global-local", 
                                                                              gapOpening = 300, substitutionMatrix = BLOSUM45)
              
              percent_identity_realignment <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
              
              print(hold_MAMP_seqs[i,])
              print(paste(hold_MAMP_seqs$Percent_Identity[i], percent_identity_realignment, sep = "_"))
              
              if (hold_MAMP_seqs$Percent_Identity[i] >= percent_identity_realignment){
                filter_list_2[[i]] <- TRUE
              }
              if (hold_MAMP_seqs$Percent_Identity[i] < percent_identity_realignment){
                hold_MAMP_seqs$MAMP_Sequence[i] <- as.character(Alignment_between_MAMP_and_Ref@subject)
                hold_MAMP_seqs$Percent_Identity[i] <- percent_identity_realignment
                filter_list_2[[i]] <- TRUE
              }
            }
          }

        }
        
        
        # flagellin - flg22
        if (hold_MAMP_seqs$MAMP_Hit[i] == "flg22_consensus"){
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flg22_consensus')
          hold_MAMP_seqs$MAMP_Sequence[i] <- str_replace_all(hold_MAMP_seqs$MAMP_Sequence[i], "-", "")
          if (hold_MAMP_seqs$Percent_Identity[i] > 30){
            filter_list_2[[i]] <- TRUE
          }
          if (hold_MAMP_seqs$Percent_Identity[i] < 30){
            print(hold_MAMP_seqs[i,1:4])
            filter_list_2[[i]] <- FALSE
          }
        }
        
        
        # flagellin - flgII-28
        if (hold_MAMP_seqs$MAMP_Hit[i] == "flgII-28"){
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flgII-28')
          hold_MAMP_seqs$MAMP_Sequence[i] <- str_replace_all(hold_MAMP_seqs$MAMP_Sequence[i], "-", "")
          if (hold_MAMP_seqs$Percent_Identity[i] > 30){
            filter_list_2[[i]] <- TRUE
          }
          if (hold_MAMP_seqs$Percent_Identity[i] < 30){
            filter_list_2[[i]] <- FALSE
          }
        }
          
        # flagellin - nlp - NPP1
        if (hold_MAMP_seqs$MAMP_Hit[i] == "nlp20_consensus"){
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'nlp20_consensus')
          hold_MAMP_seqs$MAMP_Sequence[i] <- str_replace_all(hold_MAMP_seqs$MAMP_Sequence[i], "-", "")
          if (hold_MAMP_seqs$Percent_Identity[i] > 30){
            filter_list_2[[i]] <- TRUE
          }
          if (hold_MAMP_seqs$Percent_Identity[i] < 30){
            filter_list_2[[i]] <- FALSE
          }
        }  
          
        
      }
    }
      
    
    hold_MAMP_seqs <- hold_MAMP_seqs[unlist(filter_list_2),]









# temp until I figure out what to do with nlp data

hold_copy_number <- as.data.frame(hold_MAMP_seqs %>% group_by(MAMP_Hit, File_Name, Gram, Genera) %>% summarise(n=n()))


