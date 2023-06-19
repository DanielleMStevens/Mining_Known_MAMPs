#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------





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
writeFasta(Clavibacter_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Clavibacter_CSPs.fasta")
writeFasta(Rathayibacter_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Rathayibacter_CSPs.fasta")
writeFasta(Rhodococcus_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Rhodococcus_CSPs.fasta")
writeFasta(Leifsonia_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Leifsonia_CSPs.fasta")
writeFasta(Streptomyces_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Streptomyces_CSPs.fasta")
writeFasta(Curtobacterium_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Curtobacterium_CSPs.fasta")



writeFasta(Pectobacterium_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Pectobacterium_CSPs.fasta")
writeFasta(Dickeya_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Dickeya_CSPs.fasta")
writeFasta(Erwinia_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Erwinia_CSPs.fasta")


writeFasta(Agrobacterium_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Agrobacterium_CSPs.fasta")
writeFasta(Xanthomonas_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Xanthomonas_CSPs.fasta")
writeFasta(Pseudomonas_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Pseudomonas_CSPs.fasta")
writeFasta(Ralstonia_CSPs, "./../../Analyses/Catagorizing_CSPs/CSP_types/Ralstonia_CSPs.fasta")










#--------------------------------------------------------------------------------------------------------------------------------------------------


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




