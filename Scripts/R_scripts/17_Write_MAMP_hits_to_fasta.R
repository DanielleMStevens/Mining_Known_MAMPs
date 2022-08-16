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
csp22_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")
csp22_protein_seq <- formate2fasta(csp22_protein_seq$Protein_Name, "csp22", csp22_protein_seq$Genera, csp22_protein_seq$File_Name, csp22_protein_seq$MAMP_Sequence)
writeFasta(csp22_protein_seq, "./../Protein_alignments_and_trees/csp22/csp22.fasta")

elf18_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "elf18_consensus")
elf18_protein_seq <- formate2fasta(elf18_protein_seq$Protein_Name, "elf18", elf18_protein_seq$Genera, elf18_protein_seq$File_Name, elf18_protein_seq$MAMP_Sequence)
writeFasta(elf18_protein_seq, "./../Protein_alignments_and_trees/elf18/elf18.fasta")

flg22_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flg22_consensus")
flg22_protein_seq <- formate2fasta(flg22_protein_seq$Protein_Name, "flg22", flg22_protein_seq$Genera, flg22_protein_seq$File_Name, flg22_protein_seq$MAMP_Sequence)
writeFasta(flg22_protein_seq, "./../Protein_alignments_and_trees/flg22/flg22.fasta")

flgII_28_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flgII-28")
flgII_28_protein_seq <- formate2fasta(flgII_28_protein_seq$Protein_Name, "flgII-28", flgII_28_protein_seq$Genera, flgII_28_protein_seq$File_Name, flgII_28_protein_seq$MAMP_Sequence)
writeFasta(flgII_28_protein_seq, "./../Protein_alignments_and_trees/flgII-28/flgII_28.fasta")

nlp20_protein_seq <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "nlp20_consensus")
nlp20_protein_seq <- formate2fasta(nlp20_protein_seq$Protein_Name, "nlp20_consensus", nlp20_protein_seq$Genera, nlp20_protein_seq$File_Name, nlp20_protein_seq$MAMP_Sequence)
writeFasta(nlp20_protein_seq, "./../Protein_alignments_and_trees/nlp20/nlp20.fasta")



######################################################################
# filter through blast results, filter by annotation, and put into distict fasta files
######################################################################


csp_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
extra_long_csps <- data.frame()
filC_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
EFTu_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
nlp_full_length <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

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
    
    #if csp is not in annotation (i.e, csp22 like proteins, skip over them)
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
  

  if(filtered_hold_MAMP_seqs$MAMP_Hit[i] == "nlp20_consensus"){
    
    #if csp is not in annotation (i.e, csp22 like proteins, skip over them)
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == TRUE){
      next
    }
    
    #find full length protein sequence
    if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == FALSE){
      protein_of_interest <- All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),]
      protein_of_interest <- protein_of_interest[grepl(filtered_hold_MAMP_seqs$File_Name[i], protein_of_interest$Filename),]
      
      
      nlp_full_length <- rbind(nlp_full_length, data.frame(
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
writeFasta(nlp_full_length, "./../Protein_alignments_and_trees/NLPs/nlp_full_length.fasta")


# reforte alternate "long" csps into fasta file
extra_long_csps <- formate2fasta(extra_long_csps$Protein_Name, "extra_long_csps", extra_long_csps$Genera, extra_long_csps$Filename, extra_long_csps$seq)
writeFasta(extra_long_csps, "./../Protein_alignments_and_trees/extra_long_csps.fasta")



######################################################################
# filter through fasta files and further split into genera specific protein trees
######################################################################


# seperate CSPs by genera into distinct trees
Clavibacter_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Rathayibacter_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Rhodococcus_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Leifsonia_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Streptomyces_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Curtobacterium_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

Pectobacterium_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Dickeya_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Erwinia_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

Agrobacterium_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Xanthomonas_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Pseudomonas_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Ralstonia_CSPs <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))



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
      if (protein_of_interest$Genera[1] == "Clavibacter"){
        Clavibacter_CSPs <- rbind(Clavibacter_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Xanthomonas"){
        Xanthomonas_CSPs <- rbind(Xanthomonas_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Leifsonia"){
        Leifsonia_CSPs <- rbind(Leifsonia_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Pseudomonas"){
        Pseudomonas_CSPs <- rbind(Pseudomonas_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Ralstonia"){
        Ralstonia_CSPs <- rbind(Ralstonia_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Streptomyces"){
        Streptomyces_CSPs <- rbind(Streptomyces_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Pectobacterium"){
        Pectobacterium_CSPs <- rbind(Pectobacterium_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Agrobacterium"){
        Agrobacterium_CSPs <- rbind(Agrobacterium_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Dickeya"){
        Dickeya_CSPs <- rbind(Dickeya_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Erwinia"){
        Erwinia_CSPs <- rbind(Erwinia_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Rhodococcus"){
        Rhodococcus_CSPs <- rbind(Rhodococcus_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Curtobacterium"){
        Curtobacterium_CSPs <- rbind(Curtobacterium_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Rathayibacter"){
        Rathayibacter_CSPs <- rbind(Rathayibacter_CSPs, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      
      
    }
    #print(paste(i,protein_of_interest[1,1], sep = "|"))
    
  }

  setTxtProgressBar(pb, i)
}

close(pb)


# write out data framees to fasta files 
writeFasta(Clavibacter_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Clavibacter_CSPs.fasta")
writeFasta(Rathayibacter_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Rathayibacter_CSPs.fasta")
writeFasta(Rhodococcus_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Rhodococcus_CSPs.fasta")
writeFasta(Leifsonia_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Leifsonia_CSPs.fasta")
writeFasta(Streptomyces_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Streptomyces_CSPs.fasta")
writeFasta(Curtobacterium_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Curtobacterium_CSPs.fasta")

writeFasta(Pectobacterium_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Pectobacterium_CSPs.fasta")
writeFasta(Dickeya_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Dickeya_CSPs.fasta")
writeFasta(Erwinia_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Erwinia_CSPs.fasta")

writeFasta(Agrobacterium_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Agrobacterium_CSPs.fasta")
writeFasta(Xanthomonas_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Xanthomonas_CSPs.fasta")
writeFasta(Pseudomonas_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Pseudomonas_CSPs.fasta")
writeFasta(Ralstonia_CSPs, "./../Protein_alignments_and_trees/Genera_specific_trees/Ralstonia_CSPs.fasta")


# seperate EFTu by genera for Tajima's D calculations in Figure 2 & 3
Clavibacter_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Rathayibacter_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Rhodococcus_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Leifsonia_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Streptomyces_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Curtobacterium_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

Pectobacterium_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Dickeya_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Erwinia_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))

Agrobacterium_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Xanthomonas_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Pseudomonas_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
Ralstonia_ETFu <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))


for (i in 1:nrow(filtered_hold_MAMP_seqs)){
  if(filtered_hold_MAMP_seqs$MAMP_Hit[i] == "elf18_consensus"){
    
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
      if (protein_of_interest$Genera[1] == "Clavibacter"){
        Clavibacter_ETFu <- rbind(Clavibacter_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Xanthomonas"){
        Xanthomonas_ETFu <- rbind(Xanthomonas_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Leifsonia"){
        Leifsonia_ETFu <- rbind(Leifsonia_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Pseudomonas"){
        Pseudomonas_ETFu <- rbind(Pseudomonas_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Ralstonia"){
        Ralstonia_ETFu <- rbind(Ralstonia_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Streptomyces"){
        Streptomyces_ETFu <- rbind(Streptomyces_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Pectobacterium"){
        Pectobacterium_ETFu <- rbind(Pectobacterium_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Agrobacterium"){
        Agrobacterium_ETFu <- rbind(Agrobacterium_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Dickeya"){
        Dickeya_ETFu <- rbind(Dickeya_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Erwinia"){
        Erwinia_ETFu <- rbind(Erwinia_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Rhodococcus"){
        Rhodococcus_ETFu <- rbind(Rhodococcus_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Curtobacterium"){
        Curtobacterium_ETFu <- rbind(Curtobacterium_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      if (protein_of_interest$Genera[1] == "Rathayibacter"){
        Rathayibacter_ETFu <- rbind(Rathayibacter_ETFu, data.frame(
          "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                                   protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
          "Sequence" = protein_of_interest[1,3]))
      }
      
      
    }
    #print(paste(i,protein_of_interest[1,1], sep = "|"))
    
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)


# write out data framees to fasta files 
writeFasta(Clavibacter_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Clavibacter/EFTu/Clavibacter_EFTu.fasta")
writeFasta(Rathayibacter_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Rathayibacter/EFTu/Rathayibacter_EFTu.fasta")
writeFasta(Rhodococcus_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Rhodococcus/EFTu/Rhodococcus_EFTu.fasta")
writeFasta(Leifsonia_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Leifsonia/EFTu/Leifsonia_EFTu.fasta")
writeFasta(Streptomyces_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Streptomyces/EFTu/Streptomyces_EFTu.fasta")
writeFasta(Curtobacterium_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Curtobacterium/EFTu/Curtobacterium_EFTu.fasta")

writeFasta(Pectobacterium_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Pectobacterium/EFTu/Pectobacterium_EFTu.fasta")
writeFasta(Dickeya_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Dickeya/EFTu/Dickeya_EFTu.fasta")
writeFasta(Erwinia_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Erwinia/EFTu/Erwinia_EFTu.fasta")

writeFasta(Agrobacterium_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Agrobacterium/EFTu/Agrobacterium_EFTu.fasta")
writeFasta(Xanthomonas_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Xanthomonas/EFTu/Xanthomonas_EFTu.fasta")
writeFasta(Pseudomonas_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Pseudomonas/EFTu/Pseudomonas_EFTu.fasta")
writeFasta(Ralstonia_ETFu, "./../Protein_alignments_and_trees/Genera_specific_trees/Ralstonia/EFTu/Ralstonia_EFTu.fasta")







