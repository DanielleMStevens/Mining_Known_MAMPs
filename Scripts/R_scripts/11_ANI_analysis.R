#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 12/20/2020
# Script Purpose: Processing and Plotting ANI values
# Inputs Necessary: ANI_comp_name_sheet.txt, output from fast ANI analysis
# Outputs: plots all by all comparison of genomes, organizes by clustering of ANI values, and labels based on species 
#-----------------------------------------------------------------------------------------------



########################################################################
# input files 
#########################################################################

#matrix_table <- read.table(file = file.choose()) # select the output of the ANI analysis (file named ANI_analysis)
  load_ANI_results <- list.files(path = "./../../Analyses/ANI_analysis/", pattern = ".txt_ANI_comparison")


#########################################################################
# update table variables
#########################################################################


  ANI_values_df <- data.frame("Genome1" = character(0), "Genome2" = character(0), "ANI_val" = numeric(0))
  for (i in 1:length(load_ANI_results)){
    pb <- txtProgressBar(min = 0, max = length(load_ANI_results), style = 3)
    ANI_output <- read.table(file = paste("./../../Analyses/ANI_analysis/", load_ANI_results[[i]], sep = "")) 
    colnames(ANI_output) <- c("Genome1","Genome2","ANI_val","val2","val3")
    ANI_output$Genome1 <- as.character(ANI_output$Genome1)
    ANI_output$Genome2 <- as.character(ANI_output$Genome2)
    ANI_output$val2 <- NULL
    ANI_output$val3 <- NULL
    
    
    ANI_output <- ANI_output %>% mutate(Genome1 = str_replace(Genome1, "/home/danimstevens/Documents/Mining_MAMPs/Whole_Genomes_v2/whole_genomes/", ""))
    ANI_output <- ANI_output %>% mutate(Genome2 = str_replace(Genome2, "/home/danimstevens/Documents/Mining_MAMPs/Whole_Genomes_v2/whole_genomes/", ""))
    
    
    ANI_output$Genome1 <- paste(vapply(strsplit(ANI_output$Genome1,"_"), `[`, 1, FUN.VALUE=character(1)), vapply(strsplit(ANI_output$Genome1,"_"), `[`, 2, FUN.VALUE=character(1)), sep="_")
    ANI_output$Genome2 <- paste(vapply(strsplit(ANI_output$Genome2,"_"), `[`, 1, FUN.VALUE=character(1)), vapply(strsplit(ANI_output$Genome2,"_"), `[`, 2, FUN.VALUE=character(1)), sep="_")
    
    
    ANI_values_df <- rbind(ANI_values_df, ANI_output)
    setTxtProgressBar(pb, i)
    rm(ANI_output)
  }
  close(pb)
  
  #########################################################################
  # swap out accession number for the file name so we can map on genera metasdata on it
  #########################################################################
  
  
  
  dt_datasettable <- data.table::data.table(datasettable)
  setkey(dt_datasettable, "Assembly_Accession")
  
  
  ANI_values_df$Genome1 <- dt_datasettable[.(ANI_values_df$Genome1),6]
  ANI_values_df$Genome2 <- dt_datasettable[.(ANI_values_df$Genome2),6]
  
  
  Genomes_to_check <- subset(ANI_values_df, ANI_values_df$ANI_val > 99.999)

  
  #########################################################################
  # Remove genome matches which are matches to itself
  #########################################################################
  
  filter_genomes <- list()
  for (i in 1:nrow(Genomes_to_check)){
    if (Genomes_to_check$Genome1[i] == Genomes_to_check$Genome2[i]){
     filter_genomes[[i]] <- FALSE
    }
    if (Genomes_to_check$Genome1[i] != Genomes_to_check$Genome2[i]){
      filter_genomes[[i]] <- TRUE
    }
  }
  Genomes_to_check <- Genomes_to_check[unlist(filter_genomes),]
  rm(filter_genomes)
  gc()
  
  #########################################################################
  # Remove genome matches which are matches to itself
  #########################################################################
  
  Do_MAMPs_match <- list()
  for (i in 1:nrow(Genomes_to_check)){
    Genome1_MAMPs <- hold_MAMP_seqs[hold_MAMP_seqs$File_Name %in% Genomes_to_check$Genome1[i],]
    Genome2_MAMPs <- hold_MAMP_seqs[hold_MAMP_seqs$File_Name %in% Genomes_to_check$Genome2[i],]
    
    if (all(Genome1_MAMPs$MAMP_Sequence %in% Genome2_MAMPs$MAMP_Sequence) == TRUE){
      Do_MAMPs_match[[i]] <- TRUE
    }
    
    if (all(Genome1_MAMPs$MAMP_Sequence %in% Genome2_MAMPs$MAMP_Sequence) == FALSE){
      Do_MAMPs_match[[i]] <- FALSE
    }
  }
  
  Genomes_to_check <- Genomes_to_check[unlist(Do_MAMPs_match),]
  rm(Genome1_MAMPs)
  rm(Genome2_MAMPs)

  
  # remove reciprocal comparison
  Genomes_to_check <- Genomes_to_check[!duplicated(cbind(pmin(Genomes_to_check[,1], Genomes_to_check[,2]), pmax(Genomes_to_check[,1], Genomes_to_check[,2]))),]
  
  
  # we can use this to remove all of column 1 but keep column 2, hence remove duplicates 
  Genomes_to_check <- data.frame(lapply(Genomes_to_check, unlist))
  rownames(Genomes_to_check) <- NULL
  filtered_hold_MAMP_seqs <- hold_MAMP_seqs[!hold_MAMP_seqs$File_Name %in% Genomes_to_check$Genome1,]
  
  
