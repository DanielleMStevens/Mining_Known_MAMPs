#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


##############################################
# Plot the similarity of the MAMP to the consensus seperate by MAMP and genera
##############################################


plot_points_Percent_Ident_of_MAMPs <- function(MAMP_of_interest){
  new_plot <- ggplot(MAMP_of_interest, aes(x = factor(Genera, level = name_list), y = Percent_Identity, colour = Genera, fill = Genera)) +
    see::geom_violinhalf(color = NA, trim = T, scale = "width") +
    geom_boxplot(color = "black", outlier.alpha = 0, width = 0.15) +
    scale_color_manual("Genera", values = Genera_colors) +
    scale_fill_manual("Genera", values = Genera_colors) +
    my_ggplot_theme +
    xlab("Genera\n") +
    ylab("Percent AA\nSimilarity") +
    scale_y_continuous(limits = c(0,100),
                       breaks = c(0,20,40,60,80,100)) +
    #ggtitle(name_of_MAMP) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 0),
          panel.grid.major.x = element_line(size = 0.5, color = "grey92")) 
    #coord_flip()
  
  return(new_plot)
}

plot_copy_number <- function(data_to_plot){
  
  new_plot <- ggplot(data_to_plot, aes(x = factor(Genera, level = name_list), y = n, colour = Genera, fill = Genera)) +
    see::geom_violinhalf(color = "black", trim = T, scale = "width") +
    geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.15) +
    scale_color_manual("Genera", values = Genera_colors) +
    scale_fill_manual("Genera", values = Genera_colors) +
    my_ggplot_theme +
    xlab("Genera\n") +
    ylab("Number of Protein \n Coding Epitopes") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = ceiling(max(data_to_plot$n)/1.5))) +
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 9),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(size = 0.5, color = "grey92")
    ) +
    coord_flip()
  
  return(new_plot)
  
}


# prepare and plot the data of similarity versus copy number for csp22 eptitopes
Similarity_csp22 <- plot_points_Percent_Ident_of_MAMPs(subset(filtered_hold_MAMP_seqs, MAMP_Hit == 'csp22_consensus')) 
Copy_number_csp22 <- plot_copy_number(subset(hold_copy_number, MAMP_Hit == "csp22_consensus"))
csp22_combined_plots <- Similarity_csp22 + Copy_number_csp22

ggsave(csp22_combined_plots, filename = "./../Figures/Plot_MAMP_similarity_by_genera_csp22.pdf", device = cairo_pdf, width = 4, height = 5, units = "in")
rm(csp22_combined_plots)

# prepare and plot the data of similarity versus copy number for elf18 eptitopes
Similarity_elf18 <- plot_points_Percent_Ident_of_MAMPs(subset(filtered_hold_MAMP_seqs, MAMP_Hit == 'elf18_consensus'))
Copy_number_elf18 <- plot_copy_number(subset(hold_copy_number, MAMP_Hit == "elf18_consensus"))
#hold_elf18_copy_number_data <- hold_elf18_copy_number_data[hold_elf18_copy_number_data$n != 0, ]
Similarity_elf18 + Copy_number_elf18

ggsave(elf18_combined_plots, filename = "./../Figures/Plot_MAMP_similarity_by_genera_elf18.pdf", device = cairo_pdf, width = 4.5, height = 4, units = "in")


# prepare and plot the data of similarity versus copy number for flg22 eptitopes
Similarity_flg22 <- plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'flg22_consensus')) 
hold_flg22_copy_number_data <- subset(hold_copy_number, MAMP_Hit == "flg22_consensus")
hold_flg22_copy_number_data <- hold_flg22_copy_number_data[hold_flg22_copy_number_data$n != 0, ]
Copy_number_flg22 <- plot_copy_number(hold_flg22_copy_number_data)

Similarity_flg22 + Copy_number_flg22

#ggsave(, filename = "./../Figures/Plot_MAMP_similarity_by_genera_flg22.pdf", device = cairo_pdf, width = 5.8, height = 3.8, units = "in")


# prepare and plot the data of similarity versus copy number for flgII-28 eptitopes
Similarity_flgII28 <- plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'flgII-28')) 
hold_flgII28_copy_number_data <- subset(hold_copy_number, MAMP_Hit == "flgII-28")
hold_flgII28_copy_number_data <- hold_flgII28_copy_number_data[hold_flgII28_copy_number_data$n != 0, ]
Copy_number_flgII28 <- plot_copy_number(hold_flgII28_copy_number_data)

Similarity_flgII28 + Copy_number_flgII28






