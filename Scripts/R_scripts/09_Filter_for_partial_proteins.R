#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/08/2023
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


######################################################################
# filter for parital proteins
######################################################################


hold_df_WP_partial <- All_target_by_annotation[grepl("partial", All_target_by_annotation$names),]
hold_df_WP_partial <- unique(hold_df_WP_partial$`unlist(hold_WP_tag)`)


filter_list_2 <- list()
for (i in 1:nrow(hold_MAMP_seqs)){
  if (any(grepl(hold_MAMP_seqs$Protein_Name[i], hold_df_WP_partial)) == FALSE){
    filter_list_2[[i]] <- TRUE
  }
  if (any(grepl(hold_MAMP_seqs$Protein_Name[i], hold_df_WP_partial)) == TRUE){
    filter_list_2[[i]] <- FALSE
  }
}

hold_MAMP_seqs <- hold_MAMP_seqs[unlist(filter_list_2),]
rm(filter_list_2)


######################################################################
# final mannnual filter
######################################################################

# these are hits I've manuallly assesed that slipped through other filtering parameters but are very likely not relavent
filter_list_2 <- c("WP_146034239.1","WP_027713652.1")
hold_MAMP_seqs <- hold_MAMP_seqs[!hold_MAMP_seqs$Protein_Name %in% filter_list_2,]

## these hits (WP_174040680.1) appear to have a 1 aa deletion in the center which casues issues with the alogrithum, droppping off the 
# 1 aa length of the epitope. I manually aligned the sequence and replace it below

# replace the sequence
hold_MAMP_seqs[hold_MAMP_seqs$Protein_Name %in% "WP_174040680.1",4] <- "ATGTVKFNATKGFGFIQPDDGS"

# replace the AA sim score
hold_MAMP_seqs[hold_MAMP_seqs$Protein_Name %in% "WP_174040680.1",3] <- 
  Biostrings::pid(Biostrings::pairwiseAlignment("AVGTVKWFNAEKGFGFITPDDG", "ATGTVKFNATKGFGFIQPDDGS", type = "global-local", 
                                                                substitutionMatrix = BLOSUM45), type = "PID1")



