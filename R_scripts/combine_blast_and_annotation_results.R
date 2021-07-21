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


# filter out hits that have already been found by WP?
hold_MAMP_seqs <- hold_MAMP_seqs[,c(1,2,3,6,7,8,9,10)]
filtered_list <- All_target_by_annotation[!All_target_by_annotation$Protein_Name %in% hold_MAMP_seqs$Protein_Name,]
filtered_list <- filtered_list[,c(8,11,10,9,4,5,6,7)]
rownames(filtered_list) <- NULL
colnames(filtered_list) <- colnames(hold_MAMP_seqs)
hold_MAMP_seqs <- rbind(hold_MAMP_seqs, filtered_list)  


#remove any MAMPs below 20% similarity as the hits are likely spurious and unlikely to be real MAMPs
hold_MAMP_seqs <- subset(hold_MAMP_seqs, hold_MAMP_seqs$Percent_Identity > 20)

# temp until I figure out what to do with nlp data
hold_MAMP_seqs <- subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit != "nlp20_consensus")

hold_copy_number <- as.data.frame(hold_MAMP_seqs %>% group_by(MAMP_Hit, File_Name) %>% summarise(n=n()))



######################################################################
# go through hits to make sure there are no off targets
######################################################################
