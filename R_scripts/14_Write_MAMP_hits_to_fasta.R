#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# filter through blast results, filter by annotation, and put into distict fasta files
######################################################################


# dataframe for each combo
csp22_protein_seq <- subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")
csp22_protein_seq <- formate2fasta(hold_MAMP_seqs$Protein_Name, "csp22", hold_MAMP_seqs$Genera, hold_MAMP_seqs$File_Name, hold_MAMP_seqs$MAMP_Sequence)

elf18_protein_seq <- subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "elf18_consensus")
elf18_protein_seq <- formate2fasta(hold_MAMP_seqs$Protein_Name, "elf18", hold_MAMP_seqs$Genera, hold_MAMP_seqs$File_Name, hold_MAMP_seqs$MAMP_Sequence)

flg22_protein_seq <- subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "flg22_consensus")
flg22_protein_seq <- formate2fasta(hold_MAMP_seqs$Protein_Name, "flg22", hold_MAMP_seqs$Genera, hold_MAMP_seqs$File_Name, hold_MAMP_seqs$MAMP_Sequence)

flgII_28_protein_seq <- subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "flgII-28")
flgII_28_protein_seq <- formate2fasta(hold_MAMP_seqs$Protein_Name, "flgII-28", hold_MAMP_seqs$Genera, hold_MAMP_seqs$File_Name, hold_MAMP_seqs$MAMP_Sequence)


#csp_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
#filC_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
#EFTu_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))



######################################################################
# filter through blast results, filter by annotation, and put into distict fasta files
######################################################################


writeFasta(csp22_protein_seq, "./../Protein_alignments_and_trees/csp22/csp22.fasta")
writeFasta(elf18_protein_seq, "./../Protein_alignments_and_trees/elf18/elf18.fasta")
writeFasta(flg22_protein_seq, "./../Protein_alignments_and_trees/flg22/flg22.fasta")
writeFasta(flgII_28_protein_seq, "./../Protein_alignments_and_trees/flgII-28/flgII_28.fasta")


#writeFasta(csp_full_length, "./../Protein_alignments_and_trees/cold_shock_protein/csp22_full_length.fasta")
#writeFasta(EFTu_full_length, "./../Protein_alignments_and_trees/EfTu/elf18_full_length.fasta")
#writeFasta(filC_full_length, "./../Protein_alignments_and_trees/Flagellin/flg22_full_length.fasta")




#pb <- txtProgressBar(min = 0, max = nrow(hold_MAMP_seqs), style = 3)

#for (i in 1:nrow(hold_MAMP_seqs)){
#  if(hold_MAMP_seqs$MAMP_Hit[i] == "csp22_consensus"){
#    protein_of_interest <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,]
#    if(grepl("cold", protein_of_interest[,2]) == T){
#      # put MAMP in database
#      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "MAMP_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
#                            hold_MAMP_seqs$MAMP_Sequence[i])
#      colnames(temp_df) <- colnames(csp22_protein_seq)
#      csp22_protein_seq <- rbind(csp22_protein_seq, temp_df)
#      
#      #find full length protein sequence
#      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "Full_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
#                            protein_of_interest[,3])
#      colnames(temp_df) <- colnames(csp_full_length)
#      csp_full_length <- rbind(csp_full_length, temp_df)
#    }
#  }
#  if(hold_MAMP_seqs$MAMP_Hit[i] == "flg22_consensus"){
#    protein_of_interest <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,]
#    if(grepl("Agrobacterium", protein_of_interest[,2]) == T){
#      print(protein_of_interest)
#    }
#    if(grepl("flagellin", protein_of_interest[,2]) == T){
#      
#      
#      # put MAMP in database
#      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "MAMP_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
#                            hold_MAMP_seqs$MAMP_Sequence[i])
#      colnames(temp_df) <- colnames(flg22_protein_seq)
#      flg22_protein_seq <- rbind(flg22_protein_seq, temp_df)
#      
#      #find full length protein sequence
#      pull_protein_seq <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,3]
#      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "Full_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
#                            protein_of_interest[,3])
#      colnames(temp_df) <- colnames(filC_full_length)
#      filC_full_length <- rbind(filC_full_length, temp_df)
#    }
#  }
#  if(hold_MAMP_seqs$MAMP_Hit[i] == "elf18_consensus"){
#    protein_of_interest <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,]
#    if(grepl("factor", protein_of_interest[,2]) == T){
#      
#      
#      # put MAMP in database
#      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "MAMP_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
#                            hold_MAMP_seqs$MAMP_Sequence[i])
#      colnames(temp_df) <- colnames(elf18_protein_seq)
#      elf18_protein_seq <- rbind(elf18_protein_seq, temp_df)
#      
#      #find full length protein sequence
#      pull_protein_seq <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,3]
#      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "Full_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
#                            protein_of_interest[,3])
#      colnames(temp_df) <- colnames(EFTu_full_length)
#      EFTu_full_length <- rbind(EFTu_full_length, temp_df)
#    }
#  }
  
  
#  setTxtProgressBar(pb, i)
#}

#close(pb)

