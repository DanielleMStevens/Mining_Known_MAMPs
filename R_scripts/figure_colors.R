#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 6/30/2020
# Script Purpose: Hold information of colors for all plots/figures 
# Inputs Necessary: n/a
# Outputs: n/a
#-----------------------------------------------------------------------------------------------


##############################################
# Colors for Script "Plot_PAMP_responses"
##############################################


Pathogenicity_colors <- c("Black","Grey") 
names(Pathogenicity_colors) <- c("Yes","No")


Genera_colors <- c("#ff94af","#ffd65c","#b589d6",
                   "#73bfe6","#9e9e9e","#7bc98f",
                   "#ab234c","#e3534c","#026178",
                   "#507B00")

names(Genera_colors) <- c("Leifsonia","Clavibacter","Curtobacterium",
                          "Rhodococcus","Rathayibacter","Streptomyces",
                          "Ralstonia", "Xanthomonas","Pseudomonas",
                          "Agrobacterium")


# defines a name list to consistently order genera
name_list <- c('Clavibacter','Leifsonia','Rathayibacter','Curtobacterium','Rhodococcus','Streptomyces',
               'Agrobacterium','Ralstonia','Xanthomonas','Pseudomonas')



MAMP_colors <- c("#FDBE83", "#C8A3B5", "#2F4E68")
names(MAMP_colors) <- c("csp22_consensus","elf18_consensus","flg22_consensus")


Gram_colors <- c("#204254","#B94A3E")
names(Gram_colors) <- c("+","-")