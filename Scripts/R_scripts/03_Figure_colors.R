#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 6/30/2020
# Script Purpose: Hold information of colors for all plots/figures 
# Inputs Necessary: N/A
# Outputs: N/A
#-----------------------------------------------------------------------------------------------


##############################################
# Colors for Script "Plot_PAMP_responses"
##############################################


#Pathogenicity_colors <- c("Black","Grey") 
#names(Pathogenicity_colors) <- c("Yes","No")

# defines a name list to consistently order genera
name_list <- c('Clavibacter','Leifsonia','Rathayibacter','Curtobacterium','Rhodococcus','Streptomyces',
               'Agrobacterium','Ralstonia','Xanthomonas','Pseudomonas','Pectobacterium','Dickeya','Erwinia')



# color code for genera of interest
Genera_colors <- c("#ff94af","#ffd65c","#b589d6","#73bfe6","#9e9e9e","#7bc98f",
                   "#ab234c","#e3534c","#026178","#507B00","#68217a",'#bad90a','#00bcf2')

names(Genera_colors) <- c("Leifsonia","Clavibacter","Curtobacterium","Rhodococcus","Rathayibacter","Streptomyces",
                          "Ralstonia", "Xanthomonas","Pseudomonas","Agrobacterium","Pectobacterium",'Dickeya','Erwinia')


# color code for MAMP of interest
MAMP_colors <- c("#FDBE83", "#d6a5bd", "#2F4E68", "#6b8c82", "#982e26")
names(MAMP_colors) <- c("csp22_consensus","elf18_consensus","flg22_consensus","flgII-28","nlp20_consensus")


# color code for Immunogenicity
Immunogenicity_colors <- c("#000000","#005b96","#800000", "#808080")
names(Immunogenicity_colors) <- c("Immunogenic","Slightly Immunogenic","Non-Immunogenic","Not Tested")

# color code for gram-type of interest
Gram_colors <- c("#3b7796","#e05353")
names(Gram_colors) <- c("+","-")


# color code for class of interest
class_colors <- c('#a8e6cf','#dcedc1','#ffd3b6','#ffaaa5')
names(class_colors) <- c('Actinobacteria','Gammaproteobacteria','Betaproteobacteria','Alphaproteobacteria')


# set colors for ANI plots using complex heatmap
color_species <- c("Leifsonia" = "#ff94af",
                   "Clavibacter" = "#ffd65c",
                   "Curtobacterium" = "#b589d6",
                   "Rhodococcus" = "#73bfe6",
                   "Rathayibacter" = "#9e9e9e",
                   "Streptomyces" = "#7bc98f",
                   "Ralstonia" = "#ab234c",
                   "Xanthomonas" =  "#e3534c",
                   "Pseudomonas" = "#026178",
                   "Agrobacterium" = "#507B00",
                   "Pectobacterium" = "#68217a",
                   "Dickeya" = "#bad90a",
                   "Erwinia" = "#00bcf2"
)
