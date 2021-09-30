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
  load_ANI_results <- list.files(path = "./../ANI_analysis/", pattern = ".txt_ANI_comparison")


#########################################################################
# update table variables
#########################################################################

  Genomes_to_check <- data.frame("Genome1" = character(0), "Genome2" = character(0), "ANI_val" = numeric(0))
  pb <- txtProgressBar(min = 0, max = length(load_ANI_results), style = 3)
  
  for (i in 1:length(load_ANI_results)){
    ANI_output <- read.table(file = paste("./../ANI_analysis/",load_ANI_results[[i]], sep = "")) 
    colnames(ANI_output) <- c("Genome1","Genome2","ANI_val","val2","val3")
    ANI_output$Genome1 <- as.character(ANI_output$Genome1)
    ANI_output$Genome2 <- as.character(ANI_output$Genome2)
    ANI_output$val2 <- NULL
    ANI_output$val3 <- NULL
    
    
    
    ANI_output <- filter_title(ANI_output,"/home/danimstevens/Documents/Mining_MAMPs/Whole_Genomes_v2/whole_genomes/")
    
    #########################################################################
    # swap out accession number for the file name so we can map on genera metasdata on it
    #########################################################################
    
    
    for (j in 1:nrow(ANI_output)){
      hold_genome1 <- paste(strsplit(ANI_output$Genome1[j], "_")[[1]][1],
                            strsplit(ANI_output$Genome1[j], "_")[[1]][2],
                            sep = "_")
      hold_genome2 <- paste(strsplit(ANI_output$Genome2[j], "_")[[1]][1],
                            strsplit(ANI_output$Genome2[j], "_")[[1]][2],
                            sep = "_")
      ANI_output$Genome1[j] <- datasettable[datasettable$Assembly_Accession %in% hold_genome1,6]
      ANI_output$Genome2[j] <- datasettable[datasettable$Assembly_Accession %in% hold_genome2,6]
    }
    
    ANIs_above_99 <- subset(ANI_output, ANI_output$ANI_val > 99.99)
    Genomes_to_check <- rbind(Genomes_to_check, ANIs_above_99)
    setTxtProgressBar(pb, i)
    
  }

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
  
  #########################################################################
  # Remove genome matches which are matches to itself
  #########################################################################
  
  
  

