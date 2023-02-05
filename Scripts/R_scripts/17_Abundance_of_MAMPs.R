#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------

#-----------------------Figure 1 Plots------------------------------------------------------------------------------

##############################################
# Collect abundace for each MAMP
##############################################

# Assessing abundance vs. uniqueness of csp22 eptiopes
hold_total_number <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")

# Number of combinations that occur at least 10, 100 times
number_unique <- as.data.frame(hold_total_number %>% group_by(MAMP_Sequence) %>% summarise(n=n()))

hold_abundance_data_df <- data.frame("Total number of MAMP sequences assessed" = nrow(hold_total_number),
                               "Total number of combinations of MAMP epitopes" = length(unique(hold_total_number$MAMP_Sequence)),
                               "Number of combinations that occur at least 10 times" = nrow(subset(number_unique, number_unique$n > 10)),
                               "Number of combinations that occur at least 100 times" =  nrow(subset(number_unique, number_unique$n > 100)),
                               "MAMP_Hit" = "csp22_consensus"
)


# Assessing abundance vs. uniqueness of elf18 eptiopes
hold_total_number <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "elf18_consensus")

# Number of combinations that occur at least 10, 100 times
number_unique <- as.data.frame(hold_total_number %>% group_by(MAMP_Sequence) %>% summarise(n=n()))

hold_abundance_data_df <- rbind(hold_abundance_data_df,
  data.frame("Total number of MAMP sequences assessed" = nrow(hold_total_number),
                                     "Total number of combinations of MAMP epitopes" = length(unique(hold_total_number$MAMP_Sequence)),
                                     "Number of combinations that occur at least 10 times" = nrow(subset(number_unique, number_unique$n > 10)),
                                     "Number of combinations that occur at least 100 times" =  nrow(subset(number_unique, number_unique$n > 100)),
                                     "MAMP_Hit" = "elf18_consensus")
)


# Assessing abundance vs. uniqueness of flg22 eptiopes
hold_total_number <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flg22_consensus")

# Number of combinations that occur at least 10, 100 times
number_unique <- as.data.frame(hold_total_number %>% group_by(MAMP_Sequence) %>% summarise(n=n()))

hold_abundance_data_df <- rbind(hold_abundance_data_df,
                                data.frame("Total number of MAMP sequences assessed" = nrow(hold_total_number),
                                           "Total number of combinations of MAMP epitopes" = length(unique(hold_total_number$MAMP_Sequence)),
                                           "Number of combinations that occur at least 10 times" = nrow(subset(number_unique, number_unique$n > 10)),
                                           "Number of combinations that occur at least 100 times" =  nrow(subset(number_unique, number_unique$n > 100)),
                                           "MAMP_Hit" = "flg22_consensus")
)


# Assessing abundance vs. uniqueness of flgII-28 eptiopes
hold_total_number <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flgII-28")

# Number of combinations that occur at least 10, 100 times
number_unique <- as.data.frame(hold_total_number %>% group_by(MAMP_Sequence) %>% summarise(n=n()))

hold_abundance_data_df <- rbind(hold_abundance_data_df,
                                data.frame("Total number of MAMP sequences assessed" = nrow(hold_total_number),
                                           "Total number of combinations of MAMP epitopes" = length(unique(hold_total_number$MAMP_Sequence)),
                                           "Number of combinations that occur at least 10 times" = nrow(subset(number_unique, number_unique$n > 10)),
                                           "Number of combinations that occur at least 100 times" =  nrow(subset(number_unique, number_unique$n > 100)),
                                           "MAMP_Hit" = "flgII-28")
)



##############################################
# Pool together data and plot
##############################################

hold_abundance_data_df <- reshape2::melt(hold_abundance_data_df)

hold_abundance_data_df$variable <- gsub(".", " ", hold_abundance_data_df$variable, fixed = TRUE)

hold_abundance_data_df$variable <- factor(hold_abundance_data_df$variable, levels = c("Total number of MAMP sequences assessed", 
                                                                                      "Total number of combinations of MAMP epitopes", 
                                                                                      "Number of combinations that occur at least 10 times",
                                                                                      "Number of combinations that occur at least 100 times"))


#  Abundance of each MAMP occuring
abundace_figure <- ggplot(hold_abundance_data_df) +  
  geom_linerange(aes(x = MAMP_Hit, ymin = 0, ymax = value), color = "black") +
  geom_point(aes(x = MAMP_Hit, y = value, color = variable), size = 4.5) +
  my_ggplot_theme +
  scale_y_continuous(breaks = c(0,5000,10000,15000,20000),
                     limits = c(0,24000)) +
  ylab("\nNumber of Epitopes") +
  coord_flip() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.position = "bottom", legend.direction = "vertical",
        legend.title = element_blank(), legend.text = element_text(color = "black", size = 10, family = 'Arial')) +
  scale_x_discrete(name ="MAMP", 
                   labels=c("csp22_consensus" = "csp22", 
                            "elf18_consensus" = "elf18",
                            "flg22_consensus" = "flg22",
                            "flgII-28" = "flgII-28")) +
  geom_text(aes(x = MAMP_Hit, y = value, label = value),
            position = position_dodge(width = 1), vjust = -1, color = "black", family = 'Arial', size = 3.5) 





#  Zoom in plot of abundance for common combinations
abundance_zoomed_in <- abundace_figure + ylim(0,621) + theme(axis.title.y = element_blank(),
                                                             axis.text.y = element_blank(),
                                                             axis.line = element_line(colour = "black", size = 0.4, linetype = "solid"),
                                                             legend.position = "none",
                                                             plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))


ggsave(abundace_figure, filename = "./../Figures/Plot_MAMP_abundance.pdf", device = cairo_pdf, width = 5.4, height = 4.2, units = "in")


ggsave(abundance_zoomed_in, filename = "./../Figures/Supplemental_Figure_4/Plot_MAMP_abundance_Zoom_in.pdf", device = cairo_pdf, width = 5.5, height = 2.5, units = "in")










##############################################
# Nlp20 ignore
##############################################




# Assessing abundance of csps
#hold_total_number <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "nlp20_consensus")

#"Number of combinations that occur at least 15 times"
#number_unique <- as.data.frame(hold_total_number %>% group_by(MAMP_Sequence) %>% summarise(n=n()))

#hold_abundance_data_df <- rbind(hold_abundance_data_df,
#                                data.frame("Total number of MAMP sequences assessed" = nrow(hold_total_number),
#                                           "Total number of combinations of MAMP epitopes" = length(unique(hold_total_number$MAMP_Sequence)),
#                                           "Number of combinations that occur at least 10 times" = nrow(subset(number_unique, number_unique$n > 10)),
#                                           "Number of combinations that occur at least 100 times" =  nrow(subset(number_unique, number_unique$n > 100)),
#                                           "MAMP_Hit" = "nlp20_consensus")
#)




