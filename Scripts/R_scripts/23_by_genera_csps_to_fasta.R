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





#---------------------------------------------------------- export peptides from database in genera specific manner --------------------------------------------------------

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Clavibacter",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Leifsonia",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Rathayibacter",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Rhodococcus",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Curtobacterium",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Streptomyces",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Erwinia",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Pectobacterium",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Dickeya",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Xanthomonas",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Pseudomonas",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Ralstonia",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")

write_xlsx(filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus" & filtered_hold_MAMP_seqs$Genera == "Agrobacterium",] %>% arrange(-Percent_Identity),
           "/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/csp_catagorization.xlsx")


