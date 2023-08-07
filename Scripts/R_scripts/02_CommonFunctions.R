#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/08/2023
# Script Purpose: Load custom function used in paper
# Inputs: Function Depedent
# Outputs: Function Dependent 
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

# turning AAsequences (fasta) into dataframe -> uses Biostrings functionality/structure
dss2df <- function(dss){
  return(data.frame(width = dss@ranges@width, names = names(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}


# turning AAmultiplesequence alignment into dataframe -> uses Biostrings functionality/structure
aa2df <- function(dss){
  return(data.frame(names = rownames(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}



######################################################################
#  function to turn dataframe (where one column is the name and one column is the sequence)
#   into a fasta file - version 2
######################################################################

formate2fasta <- function(WP_locus_names, sequence_type, Genera_names, File_names, sequences) {
  hold_sequences <- data.frame("Locus_Tag_Name" = character(0), "Sequence" = character(0))
  pb <- txtProgressBar(min = 0, max = length(WP_locus_names), style = 3)
  for (i in 1:length(WP_locus_names)){
    #find full length protein sequence
    temp_df <- data.frame(paste(paste(">",WP_locus_names[[i]], sep=""), sequence_type, Genera_names[[i]], File_names[[i]], i, sep = "|"),
                          sequences[[i]])
    colnames(temp_df) <- colnames(hold_sequences)
    hold_sequences <- rbind(hold_sequences, temp_df)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(hold_sequences)
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



#########################################################################
# create a filter function, which will help remove information (aka the path)
# from the name of each strain such that it can be plotted with just its strain name
#########################################################################

filter_title <- function(which_table, which_phrase){
  for (i in 1:nrow(which_table)){
    which_table[i,1] <- stringr::str_replace(which_table[i,1], which_phrase, "")
    which_table[i,2] <- stringr::str_replace(which_table[i,2], which_phrase, "")
  }
  return(which_table) 
}


#########################################################################
# create a weblogo plotting function, where take a list of imput strings of the 
# same length and use ggseqlogo to create a weblogo
#########################################################################

make_me_a_weblogo <- function(file_in){
  logo <- ggseqlogo::ggseqlogo(file_in, seq_type='aa', method = 'prob') +
    theme(panel.grid = element_blank(),
          legend.position = "none", axis.text.x = element_blank(),
          axis.text.y.left = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 12, vjust = 1),
          axis.title.y = element_text(color = "black", size = 12, vjust = 1),
          axis.ticks.x = element_blank())
  
  return(logo)
}




