
######################################################################
# filter through blast results, filter by annotation, and put into distict fasta files
######################################################################


# dataframe for each combo
csp22_protein_seq <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
csp_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

flg22_protein_seq <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
filC_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

elf18_protein_seq <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
EFTu_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

pb <- txtProgressBar(min = 0, max = nrow(hold_MAMP_seqs), style = 3)

for (i in 1:nrow(hold_MAMP_seqs)){
  if(hold_MAMP_seqs$MAMP_Hit[i] == "csp22_consensus"){
    protein_of_interest <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,]
    if(grepl("cold", protein_of_interest[,2]) == T){
      # put MAMP in database
      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "MAMP_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
                            hold_MAMP_seqs$MAMP_Sequence[i])
      colnames(temp_df) <- colnames(csp22_protein_seq)
      csp22_protein_seq <- rbind(csp22_protein_seq, temp_df)
      
      #find full length protein sequence
      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "Full_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
                            protein_of_interest[,3])
      colnames(temp_df) <- colnames(csp_full_length)
      csp_full_length <- rbind(csp_full_length, temp_df)
    }
  }
  if(hold_MAMP_seqs$MAMP_Hit[i] == "flg22_consensus"){
    protein_of_interest <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,]
    if(grepl("flagellin", protein_of_interest[,2]) == T){
      
      
      # put MAMP in database
      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "MAMP_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
                            hold_MAMP_seqs$MAMP_Sequence[i])
      colnames(temp_df) <- colnames(flg22_protein_seq)
      flg22_protein_seq <- rbind(flg22_protein_seq, temp_df)
      
      #find full length protein sequence
      pull_protein_seq <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,3]
      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "Full_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
                            protein_of_interest[,3])
      colnames(temp_df) <- colnames(filC_full_length)
      filC_full_length <- rbind(filC_full_length, temp_df)
    }
  }
  if(hold_MAMP_seqs$MAMP_Hit[i] == "elf18_consensus"){
    protein_of_interest <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,]
    if(grepl("factor", protein_of_interest[,2]) == T){
      
      
      # put MAMP in database
      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "MAMP_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
                            hold_MAMP_seqs$MAMP_Sequence[i])
      colnames(temp_df) <- colnames(elf18_protein_seq)
      elf18_protein_seq <- rbind(elf18_protein_seq, temp_df)
      
      #find full length protein sequence
      pull_protein_seq <- All_target_seqs[grepl(hold_MAMP_seqs$Protein_Name[i], All_target_seqs$names),][1,3]
      temp_df <- data.frame(paste(paste(">",hold_MAMP_seqs$Protein_Name[i], sep=""), hold_MAMP_seqs$MAMP_Hit[i], "Full_Seq", hold_MAMP_seqs$Genera[i], hold_MAMP_seqs$File_Name[i], i, sep = "|"),
                            protein_of_interest[,3])
      colnames(temp_df) <- colnames(EFTu_full_length)
      EFTu_full_length <- rbind(EFTu_full_length, temp_df)
    }
  }
  
  
  setTxtProgressBar(pb, i)
}

close(pb)


######################################################################
#  function to turn dataframe (where one column is the name and one column is the sequence)
#   into a fasta file
######################################################################

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, data[rowNum,1])
    fastaLines = c(fastaLines,data[rowNum,2])
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


writeFasta(csp_full_length, "./../Protein_alignments_and_trees/cold_shock_protein/csp22_full_length.fasta")
writeFasta(EFTu_full_length, "./../Protein_alignments_and_trees/EfTu/elf18_full_length.fasta")
writeFasta(filC_full_length, "./../Protein_alignments_and_trees/Flagellin/flg22_full_length.fasta")
