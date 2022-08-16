



#-----------------------Supplemental Figure 2 Plots------------------------------------------------------------------------------

# plots all MAMPs grouped by whether it's a gram-positive or gram-negative bacteria
n_number <- as.data.frame(filtered_hold_MAMP_seqs %>% group_by(MAMP_Hit) %>% summarise(n=n()))


nlp_MAMP_sim <- ggplot(filtered_hold_MAMP_seqs, aes(x = MAMP_Hit, y = Percent_Identity)) +
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
  #scale_y_continuous(breaks = c(0,25,50,75,100), limits = c(0,110)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  geom_text(data = n_number, 
            aes(x = MAMP_Hit, y = 110, label = n), size = 4)


ggsave(nlp_MAMP_sim, filename = "./../Figures/Supplemental_Figure_2/Plot_MAMP_Sim_with_nlp20.pdf", device = cairo_pdf, width = 6, height = 3, units = "in")

rm(nlp_MAMP_sim)


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
  

ggsave(nlp_counts_figure, filename = "./../Figures/Supplemental_Figure_2/Plot_nlp20_number_per_genome.pdf", device = cairo_pdf, width = 5.5, height = 7.2, units = "in")


# plots similairt in a intra genus manner 
nlp_similarity <- subset(filtered_hold_MAMP_seqs, MAMP_Hit == 'nlp20_consensus')
nlp_similarity <- nlp_similarity[nlp_similarity$Genera %in% c("Pectobacterium", "Dickeya", "Streptomyces"),]

ggsave(plot_points_Percent_Ident_of_MAMPs(nlp_similarity) , filename = "./../Figures/Supplemental_Figure_2/Plot_nlp20_sim_per_genera.pdf", device = cairo_pdf, width = 4.2, height = 2.2, units = "in")



