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

hold_copy_number <- as.data.frame(hold_MAMP_seqs %>% group_by(MAMP_Hit, File_Name, Gram, Genera) %>% summarise(n=n()))


######################################################################
# correct ALL-Target_seqs such that it includes peptides from filtered list
######################################################################

#filtered_list <- subset(filtered_list, filtered_list$Percent_Identity > 20)
#pb <- txtProgressBar(min = 0, max = nrow(filtered_list), style = 3)

#for (j in 1:nrow(filtered_list)){
#    
#  for (i in 1:length(load_protein_fasta_files)){
    
    # cross reference with metadata to provide genus, strain, filename info
 #   get_accession_number <- strsplit(load_protein_fasta_files[[i]], "/")[[1]][6]
#    get_accession_number <- paste(strsplit(get_accession_number,"_")[[1]][1],
 #                                 strsplit(get_accession_number,"_")[[1]][2], sep = "_")
#    get_strain_info <- subset(datasettable, Assembly_Accession == get_accession_number)
 #   
  #  if (get_strain_info$Filename == filtered_list$File_Name[j]){}
  
      # read blast resukts and protein fasta file
#      read_protein_fasta <-  dss2df(Biostrings::readAAStringSet(load_protein_fasta_files[[i]]))
      
      # filter proteins in each fasta for those with the same locus tag as the blast results
 #     grab_right_protein_seq_blast_results <- read_protein_fasta[grepl(filtered_list$Protein_Name[j], read_protein_fasta$names),]
  
#      All_target_seqs <- rbind(All_target_seqs, grab_right_protein_seq_blast_results)
  
 # }
#  setTxtProgressBar(pb, j)
  
#}




######################################################################
# go through hits to make sure there are no off targets
######################################################################
