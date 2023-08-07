#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 4/20/2020
# Script Purpose: Hold information of style of ggplots 
# Inputs Necessary: N/A
# Outputs: N/A
#-----------------------------------------------------------------------------------------------


######################################################################
# library packages need to load
######################################################################

my_ggplot_theme <- theme_classic() +
  theme(axis.title.x = element_text(size = 14, color = "black", family = "Arial", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", family = "Arial", face = "bold"),
        axis.text.x = element_text(size = 12,  color = "black", family = "Arial"),
        axis.text.y = element_text(size = 12, color = "black", family = "Arial"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        axis.line = element_line(colour = "black", size = 0.4, linetype = "solid"))
        #panel.border = element_rect(color = "black", size = 0.8))
