#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 12/11/2022
# Script Purpose: Plot MAMP dynamics
# Inputs: Data processed forom earlier scripts
# Outputs: Plots of MAMP variants from diverse bacteria
#-----------------------------------------------------------------------------------------------


#-----------------------Figure 1 Plots------------------------------------------------------------------------------


##############################################
# Plot the number of MAMP encoded genes per genome seperated by MAMP
##############################################

# copy number of MAMPs (not seperated by Gram type) ----- Figure 1C ------

MAMP_type_copy_number <- ggplot(hold_copy_number,  
                                aes(x = MAMP_Hit, y = n)) +
  my_ggplot_theme +
  geom_violin(aes(fill = MAMP_Hit), alpha = 0.9, scale = "width", trim = T) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.13) +
  ylab("Number of \nMAMP Encoding Genes") +
  scale_fill_manual("MAMP", values = MAMP_colors) +
  scale_y_continuous(limits = c(0,16),
                     breaks = c(0,2,4,6,8,10,12,14,16)) +
  scale_x_discrete(name ="\nMAMP", 
                   labels=c("csp22_consensus" = "csp22", 
                            "elf18_consensus" = "elf18",
                            "flg22_consensus" = "flg22",
                            "flgII-28" = "flgII-28",
                            "nlp20_consensus" = "nlp20")) +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", face = "bold", size = 14),
        axis.title.y = element_text(color = "black", face = "bold", size = 14),
        axis.line = element_line(colour = "black", 
                                 size = 0.4, linetype = "solid")) 


MAMP_type_copy_number



ggsave(MAMP_type_copy_number, filename = "./../Figures/Figure_1/Plot_MAMP_number_by_MAMP_type.pdf", device = cairo_pdf, width = 7, height = 3.3, units = "in")




##############################################
# Plot the similarity of the MAMP to the consensus seperate by MAMP
##############################################

# plots all MAMPs grouped by MAMP type ----- Figure 1D and Supl. Figure 1A ------

n_number <- as.data.frame(filtered_hold_MAMP_seqs %>% group_by(MAMP_Hit) %>% summarise(n=n()))


MAMP_type <- ggplot(filtered_hold_MAMP_seqs, 
                    aes(x = MAMP_Hit, y = Percent_Identity)) +
  my_ggplot_theme +
  geom_violin(aes(fill = MAMP_Hit), alpha = 0.9,scale = "width", trim = T) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.13) +
  ylab("%AA Similarity \n to consensus") +
  scale_fill_manual("MAMP", values = MAMP_colors) +
  scale_x_discrete(name ="\nMAMP", 
                   labels=c("csp22_consensus" = "csp22", 
                            "elf18_consensus" = "elf18",
                            "flg22_consensus" = "flg22",
                            "flgII-28" = "flgII-28",
                            "nlp20_consensus" = "nlp20")) +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", face = "bold", size = 14),
        axis.line = element_line(colour = "black", 
                                 size = 0.4, linetype = "solid")) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  geom_text(data = n_number, 
            aes(x = MAMP_Hit, y = 110, label = n), size = 4.5)



MAMP_type <- MAMP_type + ggplot(filtered_hold_MAMP_seqs, 
                                aes(color = MAMP_Hit, x = Percent_Identity)) + 
  scale_color_manual("MAMP", values = MAMP_colors) +
  stat_ecdf(geom = "step", size = 1.25) +
  my_ggplot_theme +
  scale_x_continuous(breaks = c(0,25,50,75,100)) +
  ylab("Cumulative probability\n") +
  xlab("%AA Similarity to consensus") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", face = "bold", size = 14),
        axis.title.y = element_text(color = "black", face = "bold", size = 14, margin = unit(c(0,0,5,0), "mm")))


# Save Plot
ggsave(MAMP_type, filename = "./../Figures/Figure_1/Figure1C_D_Plot_MAMP_similarity_by_MAMP_type.pdf", device = cairo_pdf, width = 7, height = 3.3, units = "in")




#################


# Older figure - rolling ECDF of MAMP copy number -> not necessary for final manuscript but keeping code just in case

#MAMP_type_copy_number <- MAMP_type_copy_number +  ggplot(subset(hold_copy_number, n != 0), aes(color = MAMP_Hit, x = n)) + 
#stat_ecdf(geom = "step", size = 1.25) +
#my_ggplot_theme +
#scale_color_manual("MAMP", values = MAMP_colors) +
#ylab("Cumulative probability")  +
#scale_x_continuous(limits = c(0,18), breaks = c(0,2,4,6,8,10,12,14,16,18)) +
#theme(legend.position = "none")


  