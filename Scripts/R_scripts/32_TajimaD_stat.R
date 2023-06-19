#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


######################################################################
# reimport fasta file as DNA bin to load in ape to calculate Taijma D
######################################################################

#csp22_full_length_fasta <- adegenet::fasta2DNAbin(file = file.choose())
#flg22_full_length_fasta <- adegenet::fasta2DNAbin(file = "./../Protein_alignments_and_trees/Flagellin/flg22_full_length_alignment")



EFTu_full_length_alignment_Clavibacter <- pegas::tajima.test(adegenet::fasta2DNAbin(file = "./../Protein_alignments_and_trees/Genera_specific_trees/Clavibacter/EFTu/Clavibacter_EFTu_alignment.fa"))
EFTu_full_length_alignment_Ralstonia <- pegas::tajima.test(adegenet::fasta2DNAbin(file = "./../Protein_alignments_and_trees/Genera_specific_trees/Ralstonia/EFTu/Ralstonia_EFTu_alignment.fa"))


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




########
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
#plot_tajmaD(subset(hold_taijma_D_values, MAMP_Hit == "csp22_consensus"))

