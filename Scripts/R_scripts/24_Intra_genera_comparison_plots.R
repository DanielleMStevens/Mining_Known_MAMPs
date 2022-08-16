#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------

##############################################
# calculate variation between MAMPs within a given genera
##############################################

list_of_genera <- unique(filtered_hold_MAMP_seqs$Genera)
list_of_MAMPs <- unique(filtered_hold_MAMP_seqs$MAMP_Hit)

comparing_MAMP_within_genera <- data.frame("Genera" = character(0), "MAMP_Hit" = character(0), "Similarity" = numeric(0))
  for (i in 1:length(list_of_genera)){
    for (j in 1:length(list_of_MAMPs)){
      subset_MAMPs <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$Genera == list_of_genera[i])
      subset_MAMPs <- subset(subset_MAMPs, subset_MAMPs$MAMP_Hit == list_of_MAMPs[j])
      subset_MAMPs_unique <- unique(subset_MAMPs$MAMP_Sequence)
      for (k in 1:length(subset_MAMPs_unique)){
        for (l in 2:length(subset_MAMPs_unique)){
          if (l <= length(subset_MAMPs_unique)){
            #comparing_MAMP_within_genera <- rbind(comparing_MAMP_within_genera, 
                                                  print(data.frame(
                                                    "Genera" = subset_MAMPs$Genera[1],
                                                    "MAMP_Hit" = subset_MAMPs$MAMP_Hit[1],
                                                    "Similarity" = Biostrings::pid(Biostrings::pairwiseAlignment(subset_MAMPs_unique[k], 
                                                                                                                 subset_MAMPs_unique[l]))
                                                  ))
                                                  
          }
        }
      }  
      
    }
  }



