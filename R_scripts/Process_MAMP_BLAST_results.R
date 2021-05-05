#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# set path to data
######################################################################

#setwd to where repo was cloned and maintained
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
#getSrcDirectory(function(x) {x})



##############################################
# Load colors - Set tip labels to match other figures
##############################################

#make sure to set path to the same place where the figure 
source("./package_dependencies.R")

source("./CommonFunctions.R")

######################################################################
# upload raw data 
######################################################################

# load all files
load_protein_fasta_files <- list.files(path = "./../../Whole_Genomes/protein_fasta/", pattern = ".faa")
load_MAMP_blast_result_files <- list.files(path = './../MAMP_blast_outputs/', pattern = ".faa.txt")

# filter
load_protein_fasta_files <- load_protein_fasta_files[!load_protein_fasta_files %in% load_MAMP_blast_result_files]



#fix path so can load into process script
for (i in 1:length(load_protein_fasta_files)){
  load_protein_fasta_files[i] <- paste("./../../Whole_Genomes/protein_fasta/", load_protein_fasta_files[i], sep="")
  load_MAMP_blast_result_files[i] <- paste("./../MAMP_blast_outputs/", load_MAMP_blast_result_files[i], sep="")
}


# load reference file to determine percent identity of mined MAMPs
load_reference_MAMPs_fasta <- Biostrings::readAAStringSet(filepath = "./../MAMP_database/MAMP_elicitor_list.fasta") #may come with warning messgae..ignore for now
load_reference_MAMPs_fasta <- dss2df(load_reference_MAMPs_fasta)


##############################################
# load data which allows us to change accession name to strain name and add in genus info
##############################################


tip_data <- "./../Mining_for_known_MAMPs_genome_accession_info.xlsx" #choose Mining_for_known_MAMPs_genome_accession_info.xlsx
datasettable <- readxl::read_xlsx(tip_data, col_names = T)
datasettable <- as.data.frame(datasettable[,1:7])


######################################################################
# function - add protein description to blast results from protein annotation at fasta file
######################################################################

All_target_seqs <- data.frame("width" = numeric(0), "names" = character(0), "seq" = character(0))

hold_MAMP_seqs <- data.frame("Protein_Name" = character(0), "MAMP_Hit" = character(0), "Percent_Identity" = numeric(0),
                             "E-value" = numeric(0), "MAMP_length" = numeric(0), "MAMP_Sequence" = character(0), "Genera" = character(0),
                             "Strain_Name" = character(0), "File_Name" = character(0), "Gram" = character(0))

hold_copy_number <- data.frame("Genera" = character(1007), "Strain_Name" = character(1007), 
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
                                                                      gapOpening = 100, gapExtension = 100)
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
  if (any(grepl("csp22_consensus", colnames(hold_temp))) == TRUE){
    hold_copy_number[i,3] <- hold_temp$csp22_consensus
  }
  if (any(grepl("elf18_consensus", colnames(hold_temp))) == TRUE){
    hold_copy_number[i,4] <- hold_temp$elf18_consensus
  }
  if (any(grepl("flg22_consensus", colnames(hold_temp))) == TRUE){
    hold_copy_number[i,5] <- hold_temp$flg22_consensus
  }


  hold_MAMP_seqs <- rbind(hold_MAMP_seqs, read_blast_results)
  
  setTxtProgressBar(pb, i)
}

close(pb)


hold_copy_number <- reshape2::melt(hold_copy_number)
rm(hold_temp)

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
