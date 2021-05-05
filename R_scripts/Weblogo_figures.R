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
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/")
#try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
#getSrcDirectory(function(x) {x})



##############################################
# Load colors - Set tip labels to match other figures
##############################################

source("./figure_colors.R")

source("./Theme_ggplot.R")



##############################################
# Subset to make weblogos
##############################################

# seperate each group into genera
Clav_only <- subset(hold_MAMP_seqs, Genera == "Clavibacter")
Strep_only <- subset(hold_MAMP_seqs, Genera == "Streptomyces")
Leif_only <- subset(hold_MAMP_seqs, Genera == "Leifsonia")
Pseudo_only <- subset(hold_MAMP_seqs, Genera == "Pseudomonas")
Rals_only <- subset(hold_MAMP_seqs, Genera == "Ralstonia")
Curto_only <- subset(hold_MAMP_seqs, Genera == "Curtobacterium")
Xanth_only <- subset(hold_MAMP_seqs, Genera == "Xanthomonas")
Agro_only <- subset(hold_MAMP_seqs, Genera == "Agrobacterium")
Rhodo_only <- subset(hold_MAMP_seqs, Genera == "Rhodococcus")
Rathayi_only <- subset(hold_MAMP_seqs, Genera == "Rathayibacter")



make_me_a_weblogo <- function(file_in){
  logo <- ggseqlogo::ggseqlogo(file_in, seq_type='aa', method = 'prob') +
    theme(panel.grid = element_blank(),
          legend.position = "none", axis.text.x = element_blank(),
          axis.text.y.left = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 12, vjust = 1),
          axis.title.y = element_text(color = "black", size = 12, vjust = 1),
          axis.ticks.x = element_blank())
  
  return(logo)
  
  Only_G_pos <- subset(hold_MAMP_seqs, Gram == "+")
  Only_G_neg <- subset(hold_MAMP_seqs, Gram == "-")
  
  # comparisons between Gram-negative and Gram-positives csp22
  make_me_a_weblogo(Only_G_pos[Only_G_pos$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Only_G_neg[Only_G_neg$MAMP_Hit %in% 'csp22_consensus',6])
  
  # comparisons between Gram-negative and Gram-positives elf18
  make_me_a_weblogo(Only_G_pos[Only_G_pos$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Only_G_neg[Only_G_neg$MAMP_Hit %in% 'elf18_consensus',6])
  
  # weblogos ofr csp22 and csp22-like peptides
  make_me_a_weblogo(Clav_only[Clav_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Strep_only[Strep_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Leif_only[Leif_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Curto_only[Curto_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Rhodo_only[Rhodo_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Rathayi_only[Rathayi_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    
    make_me_a_weblogo(Agro_only[Agro_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Pseudo_only[Pseudo_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Rals_only[Rals_only$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Xanth_only[Xanth_only$MAMP_Hit %in% 'csp22_consensus',6])
  
  # weblogos for elf18 and elf18-like peptides
  make_me_a_weblogo(Clav_only[Clav_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Strep_only[Strep_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Leif_only[Leif_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Curto_only[Curto_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Rhodo_only[Rhodo_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Rathayi_only[Rathayi_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    
    make_me_a_weblogo(Agro_only[Agro_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Pseudo_only[Pseudo_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Rals_only[Rals_only$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Xanth_only[Xanth_only$MAMP_Hit %in% 'elf18_consensus',6])
  
  # weblogos for elf18 and elf18-like peptides
  