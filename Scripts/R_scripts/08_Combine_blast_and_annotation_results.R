#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


######################################################################
# fix flgII-28 tag
######################################################################

# fix flgII-28 name scheme in hold_MAMP seq
for (i in 1:nrow(hold_MAMP_seqs)){
  if (hold_MAMP_seqs$MAMP_Hit[i] == "flgII.28"){
    hold_MAMP_seqs$MAMP_Hit[i] <- "flgII-28"
  }
}

for (i in 1:nrow(All_target_by_annotation)){
  if(All_target_by_annotation$MAMP_Hit[i] == "flgII.28"){
    All_target_by_annotation$MAMP_Hit[i] <- "flgII-28"
  }
}

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


# readjust hit sequence if value seems off
filter_list_2 <- list()
for (i in 1:nrow(hold_MAMP_seqs)){
      
  if (hold_MAMP_seqs$Percent_Identity[i] > 40){
    hold_MAMP_seqs$MAMP_Sequence[i] <- str_replace_all(hold_MAMP_seqs$MAMP_Sequence[i], "-", "")
    filter_list_2[[i]] <- TRUE
  }
      
  if (hold_MAMP_seqs$Percent_Identity[i] <= 40){
    hold_MAMP_seqs$MAMP_Sequence[i] <- str_replace_all(hold_MAMP_seqs$MAMP_Sequence[i], "-", "")
        
        
    # cold shock protein
    if (hold_MAMP_seqs$MAMP_Hit[i] == "csp22_consensus"){
          
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'csp22_consensus')
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
              Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                              hold_MAMP_seqs$MAMP_Sequence[i], type = "global-local", 
                                                                              gapOpening = 300, substitutionMatrix = BLOSUM45)
              if (Alignment_between_MAMP_and_Ref@score > 30){
                filter_list_2[[i]] <- TRUE
              }
              if (Alignment_between_MAMP_and_Ref@score < 30){
                filter_list_2[[i]] <- FALSE
              }
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
          
          
          
          if (hold_MAMP_seqs$Percent_Identity[i] > 30){
            filter_list_2[[i]] <- TRUE
          }
          if (hold_MAMP_seqs$Percent_Identity[i] < 30){
            #print(hold_MAMP_seqs[i,1:4])
            filter_list_2[[i]] <- FALSE
          }
        }
        
        
    # flagellin - flgII-28
    if (hold_MAMP_seqs$MAMP_Hit[i] == "flgII-28"){
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flgII-28')
          if (hold_MAMP_seqs$Percent_Identity[i] > 30){
            filter_list_2[[i]] <- TRUE
          }
          if (hold_MAMP_seqs$Percent_Identity[i] < 30){
            filter_list_2[[i]] <- FALSE
          }
        }
         
         
    # NPP1 protein
    if (hold_MAMP_seqs$MAMP_Hit[i] == "nlp20_consensus"){
          if(hold_MAMP_seqs$Protein_Name == "WP_086019202.1"){
            print(hold_MAMP_seqs[i,])
          }
          pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'nlp20_consensus')

          
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
    
#remove old variables
rm(filter_list_2)
rm(filter1)
rm(filter2)
rm(filter3)



