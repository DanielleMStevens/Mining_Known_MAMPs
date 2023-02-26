#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------

# In order to build a  phylogenetic tree for each MAMP-encoding protein, first I need to filter through all my hits and 
# make a fasta file for each protein (EF-Tu, CSPs, and FliC)

######################################################################
# filter through blast results, filter by annotation, and put into distict fasta files
######################################################################


csp_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
extra_long_csps <- data.frame()
filC_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
EFTu_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))


pb <- txtProgressBar(min = 0, max = nrow(filtered_hold_MAMP_seqs), style = 3)
for (i in 1:nrow(filtered_hold_MAMP_seqs)){
  if(filtered_hold_MAMP_seqs$MAMP_Hit[i] == "csp22_consensus"){
    
    #if csp is not in annotation (i.e, csp22 like proteins, skip over them)
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == TRUE){
      next
    }
    
    #find full length protein sequence
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == FALSE){
      protein_of_interest <- All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),]
      protein_of_interest <- protein_of_interest[grepl(filtered_hold_MAMP_seqs$File_Name[i], protein_of_interest$Filename, fixed = T),]
      
      if (is.na(protein_of_interest$width[1]) == TRUE){
        print("missing seq")
        next
      }
      if (protein_of_interest$width[1] > 180){
        extra_long_csps <- rbind(extra_long_csps, protein_of_interest)
        csp_full_length <- rbind(csp_full_length, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$width[1] < 181){
        csp_full_length <- rbind(csp_full_length, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      
      
    }
    #print(paste(i,protein_of_interest[1,1], sep = "|"))
    
  }
  
  
  if(filtered_hold_MAMP_seqs$MAMP_Hit[i] == "flg22_consensus"){
    
    #if flagella is not in annotation, skip over them
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == TRUE){
      next
    }
    
    #find full length protein sequence
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == FALSE){
      protein_of_interest <- All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),]
      protein_of_interest <- protein_of_interest[grepl(filtered_hold_MAMP_seqs$File_Name[i], protein_of_interest$Filename),]
      
      
      filC_full_length <- rbind(filC_full_length, data.frame(
        "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                 protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
        "Sequence" = protein_of_interest[1,3])
      )
    }
  }
  
  
  if(filtered_hold_MAMP_seqs$MAMP_Hit[i] == "elf18_consensus"){
    
    #if csp is not in annotation (i.e, csp22 like proteins, skip over them)
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == TRUE){
      next
    }
    
    #find full length protein sequence
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == FALSE){
      protein_of_interest <- All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),]
      protein_of_interest <- protein_of_interest[grepl(filtered_hold_MAMP_seqs$File_Name[i], protein_of_interest$Filename),]
      
      
      EFTu_full_length <- rbind(EFTu_full_length, data.frame(
        "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                 protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
        "Sequence" = protein_of_interest[1,3])
      )
    }
  }
  

  
  setTxtProgressBar(pb, i)
}

close(pb)


writeFasta(csp_full_length, "./../Protein_alignments_and_trees/cold_shock_protein/csp_full_length.fasta")
writeFasta(EFTu_full_length, "./../Protein_alignments_and_trees/EfTu/EFTu_full_length.fasta")
writeFasta(filC_full_length, "./../Protein_alignments_and_trees/Flagellin/filC_full_length.fasta")
#writeFasta(nlp_full_length, "./../Protein_alignments_and_trees/NLPs/nlp_full_length.fasta")


# reforte alternate "long" csps into fasta file
extra_long_csps <- formate2fasta(extra_long_csps$Protein_Name, "extra_long_csps", extra_long_csps$Genera, extra_long_csps$Filename, extra_long_csps$seq)
writeFasta(extra_long_csps, "./../Protein_alignments_and_trees/extra_long_csps.fasta")











######################################################################
# move MAMP hit sequences (just epitope) into fasta file
######################################################################


# dataframe for each combo
#csp22_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")
#csp22_protein_seq <- formate2fasta(csp22_protein_seq$Protein_Name, "csp22", csp22_protein_seq$Genera, csp22_protein_seq$File_Name, csp22_protein_seq$MAMP_Sequence)
#writeFasta(csp22_protein_seq, "./../Protein_alignments_and_trees/csp22/csp22.fasta")

#elf18_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "elf18_consensus")
#elf18_protein_seq <- formate2fasta(elf18_protein_seq$Protein_Name, "elf18", elf18_protein_seq$Genera, elf18_protein_seq$File_Name, elf18_protein_seq$MAMP_Sequence)
#writeFasta(elf18_protein_seq, "./../Protein_alignments_and_trees/elf18/elf18.fasta")

#flg22_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flg22_consensus")
#flg22_protein_seq <- formate2fasta(flg22_protein_seq$Protein_Name, "flg22", flg22_protein_seq$Genera, flg22_protein_seq$File_Name, flg22_protein_seq$MAMP_Sequence)
#writeFasta(flg22_protein_seq, "./../Protein_alignments_and_trees/flg22/flg22.fasta")

#flgII_28_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flgII-28")
#flgII_28_protein_seq <- formate2fasta(flgII_28_protein_seq$Protein_Name, "flgII-28", flgII_28_protein_seq$Genera, flgII_28_protein_seq$File_Name, flgII_28_protein_seq$MAMP_Sequence)
#writeFasta(flgII_28_protein_seq, "./../Protein_alignments_and_trees/flgII-28/flgII_28.fasta")


##### Ignore
#nlp20_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "nlp20_consensus")
#nlp20_protein_seq <- formate2fasta(nlp20_protein_seq$Protein_Name, "nlp20_consensus", nlp20_protein_seq$Genera, nlp20_protein_seq$File_Name, nlp20_protein_seq$MAMP_Sequence)
#writeFasta(nlp20_protein_seq, "./../Protein_alignments_and_trees/nlp20/nlp20.fasta")



