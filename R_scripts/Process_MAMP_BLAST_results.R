#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# function - add protein description to blast results from protein annotation at fasta file
######################################################################
  
  All_target_seqs <- data.frame("width" = numeric(0), "names" = character(0), "seq" = character(0))
  
  hold_MAMP_seqs <- data.frame("Protein_Name" = character(0), "MAMP_Hit" = character(0), "Percent_Identity" = numeric(0),
                               "E-value" = numeric(0), "MAMP_length" = numeric(0), "MAMP_Sequence" = character(0), "Genera" = character(0),
                               "Strain_Name" = character(0), "File_Name" = character(0), "Gram" = character(0))
  
  hold_copy_number <- data.frame("Genera" = character(1007), "Strain_Name" = character(1007), "Gram" = character(1007),
                                 "csp22_consensus" = numeric(1007), "elf18_consensus" = numeric(1007), "flg22_consensus" = numeric(1007))
    
  pb <- txtProgressBar(min = 0, max = length(load_protein_fasta_files), style = 3)
  
  for (i in 1:length(load_protein_fasta_files)){
    
    # read blast resukts and protein fasta file
    read_protein_fasta <-  dss2df(Biostrings::readAAStringSet(load_protein_fasta_files[[i]]))
    read_blast_results <- read.table(file = load_MAMP_blast_result_files[[i]], sep = '\t', header = F)
    colnames(read_blast_results) <- c("Protein_Name","MAMP_Hit","Percent_Identity","E-value","MAMP_length",
                                      "Hit_Start","Hit_End","Hit_Length","Sequence")
    
    # filter proteins in each fasta for those with the same locus tag as the blast results
    grab_right_protein_seq_blast_results <- read_protein_fasta[grepl(paste(read_blast_results$Protein_Name, collapse = "|"), 
                                                                     read_protein_fasta$names),]
    
    All_target_seqs <- rbind(All_target_seqs, grab_right_protein_seq_blast_results)
  
    # fix MAMP extension error from blast results
    for (j in 1:nrow(read_blast_results)){
      if (read_blast_results$Hit_Length[j] == read_blast_results$MAMP_length[j]){
        next
      }
      if (read_blast_results$Hit_Length[j] != read_blast_results$MAMP_length[j]){
        subset_protein_seq <- grab_right_protein_seq_blast_results[grepl(read_blast_results$Protein_Name[j], 
                                                                         grab_right_protein_seq_blast_results$names),]
        pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == read_blast_results$MAMP_Hit[j])
        Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                        subset_protein_seq$seq, type = "global-local", 
                                                                        gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM62)
        read_blast_results$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
        read_blast_results$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
      }
    }
  
    # restructure to remove unncessary info
    read_blast_results <- read_blast_results[,c(1,2,3,4,5,9)]
    
    # cross reference with metadata to provide genus, strain, filename info
    get_accession_number <- strsplit(load_protein_fasta_files[[i]], "/")[[1]][6]
    get_accession_number <- paste(strsplit(get_accession_number,"_")[[1]][1],
                                    strsplit(get_accession_number,"_")[[1]][2], sep = "_")
    get_strain_info <- subset(datasettable, Assembly_Accession == get_accession_number)
    
    
    read_blast_results <- cbind(read_blast_results, 
                                rep(get_strain_info$Genera, nrow(read_blast_results)),
                                rep(get_strain_info$Strain_Name, nrow(read_blast_results)),
                                rep(get_strain_info$Filename, nrow(read_blast_results)),
                                rep(get_strain_info$Gram, nrow(read_blast_results)))
    
    
    colnames(read_blast_results) <- c("Protein_Name","MAMP_Hit","Percent_Identity","E-value","MAMP_length",
                                      "MAMP_Sequence", "Genera", "Strain_Name", "File_Name","Gram")
    
  
    
    # components for counting copy number  
    hold_temp <- data.frame(rbind(table(read_blast_results$MAMP_Hit)))
    
    
    # count the number of copies per strain
    hold_copy_number[i,1] <- unique(read_blast_results$Genera)
    hold_copy_number[i,2] <- unique(read_blast_results$Strain_Name)
    hold_copy_number[i,3] <- unique(read_blast_results$Gram)
    if (any(grepl("csp22_consensus", colnames(hold_temp))) == TRUE){
      hold_copy_number[i,4] <- hold_temp$csp22_consensus
    }
    if (any(grepl("elf18_consensus", colnames(hold_temp))) == TRUE){
      hold_copy_number[i,5] <- hold_temp$elf18_consensus
    }
    if (any(grepl("flg22_consensus", colnames(hold_temp))) == TRUE){
      hold_copy_number[i,6] <- hold_temp$flg22_consensus
    }
  
  
    hold_MAMP_seqs <- rbind(hold_MAMP_seqs, read_blast_results)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  
  hold_copy_number <- reshape2::melt(hold_copy_number)
  rm(hold_temp)

  
  
  