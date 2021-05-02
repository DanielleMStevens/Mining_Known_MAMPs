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
source("./figure_colors.R")

source("./Theme_ggplot.R")

source("./CommonFunctions.R")


##############################################
# filter out accessions that are no longer fit withn criteria of study
##############################################

only_Elf18_MAMP <- only_Elf18_MAMP[!grepl("remove_hits", only_Elf18_MAMP$Genera),]

only_Csp22_MAMP <- only_Csp22_MAMP[!grepl("remove_hits", only_Csp22_MAMP$Genera),]
only_Flg22_MAMP <- only_Flg22_MAMP[!grepl("remove_hits", only_Flg22_MAMP$Genera),]



##############################################
# Load colors - Set tip labels to match other figures
##############################################

plot_points_Percent_Ident_of_MAMPs <- function(MAMP_of_interest, name_of_MAMP){
  new_plot <- ggplot(MAMP_of_interest, aes(x = Genera, y = Percent_Identity, colour = Genera)) +
    geom_boxplot(color = "black", outlier.alpha = 0, alpha = 0.7) +
    #geom_beeswarm(size=2,priority='random')+
    geom_jitter(width = 0.25, alpha = 0.6) +
    scale_color_manual("Genera", values = Genera_colors) +
    my_ggplot_theme +
    xlab("\n Genera") +
    ylab("Percent Identity \n") +
    scale_y_continuous(limits = c(40,100),
                       breaks = c(40,50,60,70,80,90,100)) +
    ggtitle(name_of_MAMP) +
    theme(axis.text.x = element_text(angle = 90))
  
  return(new_plot)
}



tiff("Compare_all_genomes_MAMP_peptide_to_cons.tiff", height = 6450, width = 13016.67, units='px', compression = "lzw", res = 1200)

ggarrange(plot_points_Percent_Ident_of_MAMPs(only_Elf18_MAMP, "Elf18"), 
          plot_points_Percent_Ident_of_MAMPs(only_Csp22_MAMP, "Csp22"),
          plot_points_Percent_Ident_of_MAMPs(only_Flg22_MAMP, "Flg22"), 
          nrow = 1, ncol = 3, common.legend = TRUE, legend = "none")


dev.off()