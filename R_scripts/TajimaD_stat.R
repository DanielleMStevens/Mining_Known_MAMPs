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
# reimport fasta file as DNA bin to load in ape to calculate Taijma D
######################################################################

csp22_full_length_fasta <- adegenet::fasta2DNAbin(file = file.choose())
elf18_full_length_fasta <- adegenet::fasta2DNAbin(file = "./../Protein_alignments_and_trees/EfTu/elf18_full_length_alignment")
flg22_full_length_fasta <- adegenet::fasta2DNAbin(file = "./../Protein_alignments_and_trees/Flagellin/flg22_full_length_alignment")



# defines a name list to consistently order genera
name_list <- c('Clavibacter','Leifsonia','Rathayibacter','Curtobacterium','Rhodococcus','Streptomyces',
               'Agrobacterium','Ralstonia','Xanthomonas','Pseudomonas')


hold_taijma_D_values <- data.frame("Genera" = character(0),"MAMP_Hit" = character(0), "Numeber of Sequences" = numeric(0), "D-value" = numeric(0))

for (i in 1:length(name_list)){
  hold_taijma_D_values <- rbind(hold_taijma_D_values,
                                data.frame("Genera" = name_list[i],
                                           "MAMP_Hit" = "csp22_consensus",
                                           "Number of Sequences" = nrow(csp22_full_length_fasta[grepl(name_list[i], rownames(csp22_full_length_fasta)),]),
                                           "D-value" = tajima.test(csp22_full_length_fasta[grepl(name_list[i], rownames(csp22_full_length_fasta)),])$D))
  hold_taijma_D_values <- rbind(hold_taijma_D_values,
                                data.frame("Genera" = name_list[i],
                                           "MAMP_Hit" = "elf18_consensus",
                                           "Number of Sequences" = nrow(elf18__full_length_fasta[grepl(name_list[i], rownames(elf18__full_length_fasta)),]),
                                           "D-value" = tajima.test(elf18__full_length_fasta[grepl(name_list[i], rownames(elf18__full_length_fasta)),])$D))
  #hold_taijma_D_values <- rbind(hold_taijma_D_values,
  #                              data.frame("Genera" = name_list[i],
  #                                         "MAMP_Hit" = "flg22_consensus",
  #                                         "Number of Sequences" = nrow(flg22_full_length_fasta[grepl(name_list[i], rownames(flg22_full_length_fasta)),]),
  #                                        "D-value" = tajima.test(flg22_full_length_fasta[grepl(name_list[i], rownames(flg22_full_length_fasta)),])$D))
}


