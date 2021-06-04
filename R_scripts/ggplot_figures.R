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
# Plot the similarity of the MAMP to the consensus seperate by Gram-type, MAMP, or both
##############################################

## NOTE! Nlp20 data is removed for now (at least for )

  # plots all MAMPs grouped by whether it's a gram-positive or gram-negative bacteria
  Gram_type <- ggplot(subset(hold_MAMP_seqs, MAMP_Hit != "nlp20_consensus"), aes(x = Gram, y = Percent_Identity)) +
  my_ggplot_theme +
  geom_violin(fill = "grey20", alpha = 0.3, scale = "width", trim = F) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
  ylab("Percent AA Similarity") +
  theme(axis.text.x = element_text(color = "black", size = 16),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title = element_text(color = "black", face = "bold", size = 14),
          axis.line = element_line(colour = "black", 
                                   size = 0.4, linetype = "solid")) +
    ylim(0,100)
  
  ggsave(Gram_type, filename = "./Figures/Plot_MAMP_similarity_by_Gram_type.pdf", 
         device = cairo_pdf, width = 3, height = 3, units = "in")


  
##############################################
# Plot the similarity of the MAMP to the consensus seperate by Gram-type, MAMP, or both
##############################################
  
  # plots all MAMPs grouped by MAMP type
  location_color <- c("csp22_consensus" = "#00AFBB", "elf18_consensus" = "#00AFBB",
                      "flg22_consensus" = "#E7B800")
  
  

  MAMP_type <- ggplot(subset(hold_MAMP_seqs, MAMP_Hit != "nlp20_consensus"), 
                      aes(x = MAMP_Hit, y = Percent_Identity, fill = MAMP_Hit)) +
  my_ggplot_theme +
  geom_violin(scale = "width", trim = F) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
  ylab("Percent AA Similarity") +
  scale_x_discrete(name ="MAMP", 
                   labels=c("csp22_consensus" = "csp22", 
                            "elf18_consensus" = "elf18",
                            "flg22_consensus" = "flg22",
                            'nlp20_consensus' = "nlp20")) +
  scale_color_manual(values = location_color) +
  theme(axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", face = "bold", size = 14),
        axis.line = element_line(colour = "black", 
                                 size = 0.4, linetype = "solid")) +
  ylim(0,100)

  MAMP_type_copy_number <- ggplot(subset(hold_copy_number, value != 0), aes(x = variable, y = value)) +
    my_ggplot_theme +
    geom_violin(fill = "grey20", alpha = 0.3, scale = "width", trim = F) +
    geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
    ylab("Number of MAMP Encoding\nProtein Sequences") +
    scale_y_continuous(limits = c(0,14),
                       breaks = c(0,2,4,6,8,10,12,14)) +
    scale_x_discrete(name ="MAMP", 
                     labels=c("csp22_consensus" = "csp22", 
                              "elf18_consensus" = "elf18",
                              "flg22_consensus" = "flg22",
                              'nlp20_consensus' = "nlp20")) +
    theme(axis.text.x = element_text(color = "black", size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title = element_text(color = "black", face = "bold", size = 14),
          axis.line = element_line(colour = "black", 
                                   size = 0.4, linetype = "solid")) 
  
  #MAMP_type / MAMP_type_copy_number
  
  
  #ggsave(, filename = "./Figures/Plot_MAMP_similarity_by_MAMP.pdf", device = cairo_pdf, width = 3, height = 3, units = "in")


  # plots all MAMPs grouped by MAMP type and further seperated by Gram type
  
  MAMP_vs_Gram <- ggplot(subset(hold_MAMP_seqs, MAMP_Hit != "nlp20_consensus"), aes(x = MAMP_Hit, y = Percent_Identity, fill = Gram)) +
    theme_classic() +
    geom_violin(position = dodge, scale = "width") +
    geom_boxplot(aes(group = interaction(Gram, MAMP_Hit)), position = dodge, fill = "white", color = "black", outlier.alpha = 0, width = 0.1) +
    ylab("Percent AA Similarity") +
    scale_x_discrete(name ="MAMP", 
                     labels=c("csp22_consensus" = "csp22", 
                              "elf18_consensus" = "elf18",
                              "flg22_consensus" = "flg22",
                              'nlp20_consensus' = "nlp20")) +
    theme(axis.text.x = element_text(color = "black", size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title = element_text(color = "black", face = "bold", size = 14),
          axis.line = element_line(colour = "black", 
                                   size = 0.4, linetype = "solid")) +
    ylim(0,100)
  
  
  MAMP_vs_Gram_copy_number <- ggplot(subset(hold_copy_number, value != 0), aes(x = variable, y = value, fill = Gram)) +
    my_ggplot_theme +
    geom_violin(position = dodge, scale = "width") +
    geom_boxplot(aes(group = interaction(Gram, variable)), position = dodge, fill = "white", color = "black", outlier.alpha = 0, width = 0.1) +
    ylab("Number of MAMP Encoding\nProtein Sequences") +
    scale_y_continuous(limits = c(0,14),
                       breaks = c(0,2,4,6,8,10,12,14)) +
    scale_x_discrete(name ="MAMP", 
                     labels=c("csp22_consensus" = "csp22", 
                              "elf18_consensus" = "elf18",
                              "flg22_consensus" = "flg22",
                              'nlp20_consensus' = "nlp20")) +
    theme(axis.text.x = element_text(color = "black", size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title = element_text(color = "black", face = "bold", size = 14),
          axis.line = element_line(colour = "black", 
                                   size = 0.4, linetype = "solid")) 
  
  (MAMP_type | MAMP_vs_Gram) + plot_layout(guides = 'collect')
  
  #/ (MAMP_type_copy_number | MAMP_vs_Gram_copy_number)
  
  #ggsave(MAMP_vs_Gram, filename = "./Figures/Plot_MAMP_similarity_by_MAMP_and_Gram_type.pdf", device = cairo_pdf, width = 4.5, height = 3, units = "in")
  



##############################################
# Plot the similarity of the MAMP to the consensus seperate by MAMP and genera
##############################################

  # defines a name list to consistently order genera
  name_list <- c('Clavibacter','Leifsonia','Rathayibacter','Curtobacterium','Rhodococcus','Streptomyces',
                 'Agrobacterium','Ralstonia','Xanthomonas','Pseudomonas')
  
  
  plot_points_Percent_Ident_of_MAMPs <- function(MAMP_of_interest){
    new_plot <- ggplot(MAMP_of_interest, aes(x = factor(Genera, level = name_list), y = Percent_Identity, colour = Genera, fill = Genera)) +
      geom_violinhalf(color = "black", trim = T, scale = "width") +
      geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.15) +
      scale_color_manual("Genera", values = Genera_colors) +
      scale_fill_manual("Genera", values = Genera_colors) +
      my_ggplot_theme +
      xlab("Genera\n") +
      ylab("\nPercent AA Similarity") +
      scale_y_continuous(limits = c(0,100),
                         breaks = c(0,20,40,60,80,100)) +
      #ggtitle(name_of_MAMP) +
      theme(legend.position = "none",panel.grid.major.x = element_line(size = 0.5, color = "grey92")) +
      coord_flip()
    
    return(new_plot)
  }
  
  plot_copy_number <- function(data_to_plot){
    
    new_plot <- ggplot(data_to_plot, aes(x = factor(Genera, level = name_list), y = value, colour = Genera, fill = Genera)) +
      geom_violinhalf(color = "black", trim = T, scale = "width") +
      geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.15) +
      scale_color_manual("Genera", values = Genera_colors) +
      scale_fill_manual("Genera", values = Genera_colors) +
      my_ggplot_theme +
      xlab("Genera\n") +
      ylab("\nNumber of Protein \n Coding Epitopes") +
      scale_y_continuous(breaks= pretty_breaks()) +
      theme(legend.position = "none", 
            axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line(size = 0.5, color = "grey92")
             ) +
      coord_flip()
    
    return(new_plot)
    
  }
  
  plot_tajmaD <- function(data_to_plot){
    new_plot <- ggplot(data_to_plot, aes(x = factor(Genera, level = name_list), y = D.value, colour = Genera, fill = Genera)) +
      geom_point(size = 4, shape = 21, colour = "black") +
      scale_color_manual("Genera", values = Genera_colors) +
      scale_fill_manual("Genera", values = Genera_colors) +
      scale_y_continuous(limits = c(-4,4),
                         breaks = c(-4,-3,-2,-1,0,1,2,3,4),
                         labels = c(-4,-3,-2,-1,0,1,2,3,4)) +
      ylab("\nD-value") +
      xlab("Genera\n") +
      my_ggplot_theme +
      theme(legend.position = "none", 
            axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line(size = 0.5, color = "grey92")) +
      coord_flip() 
   
    return(new_plot) 
  }
  
  
  
  # prepare and plot the data of similarity versus copy number for csp22 eptitopes
  Similarity_csp22 <- plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'csp22_consensus')) 
  Copy_number_csp22 <- plot_copy_number(subset(hold_copy_number, variable == "csp22_consensus"))
  csp22_combined_plots <- Similarity_csp22 + Copy_number_csp22 + plot_tajmaD(subset(hold_taijma_D_values, MAMP_Hit == "csp22_consensus"))
  
  ggsave(csp22_combined_plots, filename = "./../Figures/Plot_MAMP_similarity_by_genera_csp22.pdf", 
         device = cairo_pdf, width = 8, height = 5, units = "in")
  
  
  # prepare and plot the data of similarity versus copy number for elf18 eptitopes
  Similarity_elf18 <- plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'elf18_consensus'))
  hold_elf18_copy_number_data <- subset(hold_copy_number, variable == "elf18_consensus")
  hold_elf18_copy_number_data <- hold_elf18_copy_number_data[hold_elf18_copy_number_data$value != 0, ]
  Copy_number_elf18 <- plot_copy_number(hold_elf18_copy_number_data)
  elf18_combined_plots <- Similarity_elf18 + Copy_number_elf18 + plot_tajmaD(subset(hold_taijma_D_values, MAMP_Hit == "elf18_consensus"))
    
  ggsave(elf18_combined_plots, filename = "./../Figures/Plot_MAMP_similarity_by_genera_elf18.pdf", 
         device = cairo_pdf, width = 8, height = 5, units = "in")
  
  
  # prepare and plot the data of similarity versus copy number for elf18 eptitopes
  Similarity_flg22 <- plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'flg22_consensus')) 
  hold_flg22_copy_number_data <- subset(hold_copy_number, variable == "flg22_consensus")
  hold_flg22_copy_number_data <- hold_flg22_copy_number_data[hold_flg22_copy_number_data$value != 0, ]
  Copy_number_flg22 <- plot_copy_number(hold_flg22_copy_number_data)
    
  
  ggsave(Similarity_flg22 + Copy_number_flg22, filename = "./../Figures/Plot_MAMP_similarity_by_genera_flg22.pdf", device = cairo_pdf, width = 5.8, height = 3.8, units = "in")
  
  #plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'nlp20_consensus'), "nlp20")
  








