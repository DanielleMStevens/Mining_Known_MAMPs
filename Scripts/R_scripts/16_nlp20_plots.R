#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


#-----------------------Supplemental Figure 1A Plots------------------------------------------------------------------------------


# plots all MAMPs grouped by whether it's a gram-positive or gram-negative bacteria
n_number <- as.data.frame(filtered_hold_MAMP_seqs %>% group_by(MAMP_Hit) %>% summarise(n=n()))


nlp_MAMP_sim <- ggplot(filtered_hold_MAMP_seqs, aes(x = MAMP_Hit, y = Percent_Identity)) +
  my_ggplot_theme +
  geom_violin(aes(fill = MAMP_Hit), alpha = 0.9,scale = "width", trim = T) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.15) +
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
  #scale_y_continuous(breaks = c(0,25,50,75,100), limits = c(0,110)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  geom_text(data = n_number, 
            aes(x = MAMP_Hit, y = 110, label = n), size = 4)




#-----------------------Supplemental Figure 1B Plots------------------------------------------------------------------------------

# plots percent of nlp hits per genus of bacteria
nlp_colors <- c("#B4B4B4", "#494848")
names(nlp_colors) <- c("number_of_genomes", "genomes_with_nlp20")

nlp_counts <- as.data.frame(subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "nlp20_consensus") %>% group_by(Genera) %>% summarise('genomes_with_nlp20'= n()))
number_of_genomes <- as.data.frame(set_colors_df[!set_colors_df$Filename %in% Genomes_to_check$Genome1,] %>% group_by(Genera) %>% summarise('number_of_genomes' = n()))

nlp_counts <- merge(number_of_genomes, nlp_counts, all=TRUE)
nlp_counts[is.na(nlp_counts)] <- 0

nlp_counts <- reshape2::melt(nlp_counts)

nlp_counts_figure <- ggplot(nlp_counts, aes(x = Genera, y = value)) +   
  geom_bar(aes(fill = variable), position = "dodge", stat="identity", color = "black", size = 0.3) +
  my_ggplot_theme +
  coord_flip() +
  scale_fill_manual("", label = c("Number of Genomes", "Number of nlp20 Hits"),  values = nlp_colors) +
  ylab("") +
  xlab("Genera\n") +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14, face = 'bold', family = 'Arial'),
        legend.text = element_text(color = "black", size = 12, family = 'Arial'),
        legend.title = element_text(color = "black", size = 12, face = 'bold', family = 'Arial'),
        legend.position = 'bottom',
        legend.direction = "vertical",
        axis.line = element_line(colour = "black", 
                                 size = 0.4, linetype = "solid")) +
  scale_y_continuous(breaks = c(0,250,500,750,1000,1250,1500),
                     limits = c(0,1600)) +
  geom_text(aes(x = Genera, y = value, label = value, group = variable),
    position = position_dodge(width = 1), size = 4, hjust = -0.3, color = "black", family = 'Arial')
  





#-----------------------Supplemental Figure 1C Plots------------------------------------------------------------------------------


# plots similairt in a intra genus manner 
nlp_similarity <- subset(filtered_hold_MAMP_seqs, MAMP_Hit == 'nlp20_consensus')
nlp_similarity <- nlp_similarity[nlp_similarity$Genera %in% c("Pectobacterium", "Dickeya", "Streptomyces"),]

plot_points_Percent_Ident_of_MAMPs(nlp_similarity) + theme(axis.text.x = element_text(size = 11),
                                                           axis.text.y = element_text(size = 11),
                                                           axis.title.x = element_text(size = 12),
                                                           axis.title.y = element_text(size = 12))









############## OLD CODE - IGNORE FOr NOW ############################


if(filtered_hold_MAMP_seqs$MAMP_Hit[i] == "nlp20_consensus"){
  
  #if csp is not in annotation (i.e, csp22 like proteins, skip over them)
  if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == TRUE){
    next
  }
  
  #find full length protein sequence
  if (any(is.na(All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),][1,])) == FALSE){
    protein_of_interest <- All_target_by_annotation[grepl(filtered_hold_MAMP_seqs$Protein_Name[i], All_target_by_annotation$names),]
    protein_of_interest <- protein_of_interest[grepl(filtered_hold_MAMP_seqs$File_Name[i], protein_of_interest$Filename),]
    
    
    nlp_full_length <- rbind(nlp_full_length, data.frame(
      "Locus_Tag_Name" = paste(paste(">",protein_of_interest$Protein_Name[1], sep=""), protein_of_interest$MAMP_Hit[1], "Full_Seq", 
                               protein_of_interest$Genera[1], protein_of_interest$Filename[1], i, sep = "|"),
      "Sequence" = protein_of_interest[1,3])
    )
  }
}




