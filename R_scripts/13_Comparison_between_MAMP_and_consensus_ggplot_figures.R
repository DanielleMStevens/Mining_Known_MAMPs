#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


##############################################
# Plot the similarity of the MAMP to the consensus seperate by Gram-type, MAMP, or both
##############################################

## NOTE! Nlp20 data is removed for now (at least for )

  # plots all MAMPs grouped by whether it's a gram-positive or gram-negative bacteria
  n_number <- as.data.frame(hold_MAMP_seqs %>% group_by(Gram) %>% summarise(n=n()))

  
  Gram_type <- ggplot(data = hold_MAMP_seqs, aes(x = Gram, y = Percent_Identity)) +
  my_ggplot_theme +
  geom_violin(aes(fill = Gram), alpha = 0.9, scale = "width", trim = F) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
  scale_fill_manual("Gram Type", values = Gram_colors) +
  ylab("Percent AA Similarity") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", size = 16),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title = element_text(color = "black", face = "bold", size = 14),
          axis.line = element_line(colour = "black", 
                                   size = 0.4, linetype = "solid")) +
    scale_y_continuous(breaks = c(0,25,50,75,100),
                       limits = c(0,110)) +
    geom_text(data = n_number, aes(x = Gram, y = 110, label = n), size = 5)
  
  
  
  Gram_type + 
    ggplot(hold_MAMP_seqs, aes(color = Gram, x = Percent_Identity)) + 
    stat_ecdf(geom = "step", size = 1.5) +
    my_ggplot_theme +
    ylab("Cumulative probability") +
    scale_color_manual("Gram Type", values = Gram_colors) +
    theme(legend.position = "none")
  
  
  
  # save plot
  #ggsave(Gram_type, filename = "./../Figures/Plot_MAMP_similarity_by_Gram_type.pdf", device = cairo_pdf, width = 3, height = 3, units = "in")

  

  
##############################################
# Plot the similarity of the MAMP to the consensus seperate by MAMP
##############################################
  
  # plots all MAMPs grouped by MAMP type
 
  n_number <- as.data.frame(hold_MAMP_seqs %>% group_by(MAMP_Hit) %>% summarise(n=n()))

  
  MAMP_type <- ggplot(hold_MAMP_seqs, 
                      aes(x = MAMP_Hit, y = Percent_Identity)) +
  my_ggplot_theme +
  geom_violin(aes(fill = MAMP_Hit), alpha = 0.9,scale = "width", trim = T) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
  ylab("Percent AA Similarity") +
  scale_fill_manual("MAMP", values = MAMP_colors) +
  scale_x_discrete(name ="MAMP", 
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
    scale_y_continuous(breaks = c(0,25,50,75,100),
                       limits = c(0,110)) +
    geom_text(data = n_number, 
              aes(x = MAMP_Hit, y = 110, label = n), size = 5)
  
 
  MAMP_type <- MAMP_type + ggplot(hold_MAMP_seqs, aes(color = MAMP_Hit, x = Percent_Identity)) + 
    #stat_qq() +
    scale_color_manual("MAMP", values = MAMP_colors) +
    stat_ecdf(geom = "step", size = 1.25) +
    my_ggplot_theme +
    scale_x_continuous(breaks = c(0,25,50,75,100)) +
    ylab("Cumulative probability") +
    xlab("Percent AA Similarity") +
    theme(legend.position = "none", axis.title.y = element_text(margin = unit(c(0,0,5,0), "mm"))) 
  
  
  #ggsave(MAMP_type, filename = "./../Figures/Plot_MAMP_similarity_by_MAMP_type.pdf", device = cairo_pdf, width = 5, height = 4, units = "in")
  
  
  # copy number of MAMPs (not seperated by Gram type)
  MAMP_type_copy_number <- ggplot(subset(hold_copy_number, n != 0), aes(x = MAMP_Hit, y = n)) +
    my_ggplot_theme +
    geom_violin(aes(fill = MAMP_Hit), alpha = 0.9, scale = "width", trim = T) +
    geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
    ylab("Number of MAMP Encoding\nProtein Sequences") +
    scale_fill_manual("MAMP", values = MAMP_colors) +
    scale_y_continuous(limits = c(0,18),
                       breaks = c(0,2,4,6,8,10,12,14,16,18)) +
    scale_x_discrete(name ="MAMP", 
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
                                   size = 0.4, linetype = "solid")) 
    #geom_text(data = n_number, aes(x = MAMP_Hit, y = 18, label = n), size = 5)
  
  
    MAMP_type_copy_number <- MAMP_type_copy_number +  ggplot(subset(hold_copy_number, n != 0), aes(color = MAMP_Hit, x = n)) + 
      #stat_qq() +
      stat_ecdf(geom = "step", size = 1.25) +
      my_ggplot_theme +
      scale_color_manual("MAMP", values = MAMP_colors) +
      ylab("Cumulative probability")  +

      scale_x_continuous(limits = c(0,18),
                         breaks = c(0,2,4,6,8,10,12,14,16,18)) +
      theme(legend.position = "none")
    
  

    #ggsave(MAMP_type_copy_number, filename = "./../Figures/Plot_MAMP_number_by_MAMP_type.pdf", device = cairo_pdf, width = 5, height = 4, units = "in")
  
    
    ##############################################
    # Plot the similarity of the MAMP to the consensus seperate by Gram-type and MAMP
    ##############################################
  
    
  # plots all MAMPs grouped by MAMP type and further seperated by Gram type
  dodge <- position_dodge(width = 1)
  
  n_number <- as.data.frame(hold_MAMP_seqs %>% group_by(MAMP_Hit, Gram) %>% summarise(n=n()))
  
  
  MAMP_vs_Gram <- ggplot(hold_MAMP_seqs, aes(x = MAMP_Hit, y = Percent_Identity, fill = Gram)) +
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
    scale_y_continuous(breaks = c(0,25,50,75,100),
                       limits = c(0,110)) +
    geom_text(data = n_number, 
              aes(group = interaction(Gram, MAMP_Hit), y = 110, label = n), position = dodge, size = 5)
  
  
  #ggsave(MAMP_vs_Gram, filename = "./../Figures/Plot_MAMP_similarity_by_MAMP and Gram.pdf", device = cairo_pdf, width = 6, height = 3.5, units = "in")
  
  
  
  
  MAMP_vs_Gram_copy_number <- ggplot(subset(hold_copy_number, n != 0), aes(x = MAMP_Hit, y = n, fill = Gram)) +
    my_ggplot_theme +
    geom_violin(position = dodge, scale = "width") +
    geom_boxplot(aes(group = interaction(Gram, MAMP_Hit)), position = dodge, fill = "white", color = "black", outlier.alpha = 0, width = 0.1) +
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
  
  #ggsave(MAMP_vs_Gram_copy_number, filename = "./../Figures/Plot_MAMP_number_by_MAMP and Gram.pdf", 
  #       device = cairo_pdf, width = 6, height = 3.2, units = "in")
  
  #(MAMP_type | MAMP_vs_Gram) + plot_layout(guides = 'collect')
  
  #/ (MAMP_type_copy_number | MAMP_vs_Gram_copy_number)
  
  #ggsave(MAMP_vs_Gram, filename = "./Figures/Plot_MAMP_similarity_by_MAMP_and_Gram_type.pdf", device = cairo_pdf, width = 4.5, height = 3, units = "in")
  

  
  ##############################################
  # Filter and Plot the similarity of the MAMP to the consensus seperate by bacterial class
  ##############################################
  
  class_type <- list()
  for (i in 1:nrow(hold_MAMP_seqs)){
    class_type[[i]] <- datasettable[datasettable$Filename %in% hold_MAMP_seqs$File_Name[i],7]
    print(class_type[[i]])
  }

  class_type <- unlist(class_type)
  
  
  ################################################
  
  n_number <- as.data.frame(subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "csp22_consensus") %>% group_by(File_Name) %>% summarise(n=n()))
  
  
  ggplot(subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "csp22_consensus"), 
                      aes(x = File_Name, y = Percent_Identity)) +
    my_ggplot_theme +
    geom_jitter(aes(color = Genera), size = 0.5) +
    #geom_violin(aes(fill = MAMP_Hit), alpha = 0.9,scale = "width", trim = T) +
    #geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
    ylab("Percent AA Similarity") +
    scale_color_manual("Genera", values = Genera_colors) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title = element_text(color = "black", face = "bold", size = 14),
          axis.line = element_line(colour = "black", 
                                   size = 0.4, linetype = "solid")) +
    scale_y_continuous(breaks = c(0,25,50,75,100),
                       limits = c(0,110)) 
  

    
    
    n_number <- as.data.frame(subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "flg22_consensus") %>% group_by(File_Name) %>% summarise(n=n()))
    
    
    ggplot(subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "flg22_consensus"), 
           aes(x = File_Name, y = Percent_Identity)) +
      my_ggplot_theme +
      geom_jitter(aes(color = Genera), size = 0.5) +
      #geom_violin(aes(fill = MAMP_Hit), alpha = 0.9,scale = "width", trim = T) +
      #geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
      ylab("Percent AA Similarity") +
      scale_color_manual("Genera", values = Genera_colors) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black", size = 14),
            axis.title = element_text(color = "black", face = "bold", size = 14),
            axis.line = element_line(colour = "black", 
                                     size = 0.4, linetype = "solid")) +
      scale_y_continuous(breaks = c(0,25,50,75,100),
                         limits = c(0,110)) 
    
    geom_text(data = n_number, 
              aes(x = File_Name, y = 110, label = n), size = 5)
    
    
    ggplot(subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "flgII-28"), 
           aes(x = File_Name, y = Percent_Identity)) +
      my_ggplot_theme +
      geom_jitter(aes(color = Genera), size = 0.5) +
      #geom_violin(aes(fill = MAMP_Hit), alpha = 0.9,scale = "width", trim = T) +
      #geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
      ylab("Percent AA Similarity") +
      scale_color_manual("Genera", values = Genera_colors) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black", size = 14),
            axis.title = element_text(color = "black", face = "bold", size = 14),
            axis.line = element_line(colour = "black", 
                                     size = 0.4, linetype = "solid")) +
      scale_y_continuous(breaks = c(0,25,50,75,100),
                         limits = c(0,110)) 
  