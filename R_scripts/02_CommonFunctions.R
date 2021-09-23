#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


######################################################################
#function to adjust size 
######################################################################
# adjust size in plot plane 
plot_adjust_size <- function(desired_dpi){
  orginal_width <- dev.size('px')[1]
  orginal_height <- dev.size('px')[2]
  
  adjusted_width <- (orginal_width*desired_dpi)/72
  adjusted_height <- (orginal_height*desired_dpi)/72
  return(list(paste0("width", adjusted_width, sep = " "), paste0("height", adjusted_height), sep = " "))
}


#########################################################
# funciton - convert DNAstringset attribute to dataframe
#########################################################

# turning AAsequences (fasta) into dataframe
dss2df <- function(dss){
  return(data.frame(width = dss@ranges@width, names = names(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}


# turning AAmultiplesequence alignment into dataframe
aa2df <- function(dss){
  return(data.frame(names = rownames(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}


######################################################################
#  function to turn dataframe (where one column is the name and one column is the sequence)
#   into a fasta file
######################################################################

writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, data[rowNum,1])
    fastaLines = c(fastaLines,data[rowNum,2])
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}



######################################################################
#  function to turn dataframe (where one column is the name and one column is the sequence)
#   into a fasta file - version 2
######################################################################


formate2fasta <- function(WP_locus_names, Genera_names, File_names, sequences) {
  hold_sequences <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
  for (i in 1:length(WP_locus_names)){
    #find full length protein sequence
    temp_df <- data.frame(paste(paste(">",WP_locus_names[[i]], sep=""), Genera_names[[i]], File_names[[i]], sep = "|"),
                          sequences[[i]])
    colnames(temp_df) <- colnames(hold_sequences)
    hold_sequences <- rbind(hold_sequences, temp_df)
  }
  return(hold_sequences)
}

