
#-----------------------Supplemental Figure 1 Plots------------------------------------------------------------------------------




##############################################
# Plot the similarity of the MAMP to the consensus seperate by Gram-type and MAMP
##############################################


# plots all MAMPs grouped by MAMP type and further seperated by Gram type
dodge <- position_dodge(width = 1)

n_number <- as.data.frame( subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit != "nlp20_consensus") %>% group_by(MAMP_Hit, Gram) %>% summarise(n=n()))

MAMP_vs_Gram <- ggplot(subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit != "nlp20_consensus"), 
                       aes(x = MAMP_Hit, y = Percent_Identity, fill = Gram)) +
  theme_classic() +
  geom_violin(position = dodge, scale = "width") +
  geom_boxplot(aes(group = interaction(Gram, MAMP_Hit)), position = dodge, fill = "white", color = "black", outlier.alpha = 0, width = 0.1) +
  ylab("Percent AA Similarity") +
  scale_x_discrete(name ="MAMP", 
                   labels=c("csp22_consensus" = "csp22", 
                            "elf18_consensus" = "elf18",
                            "flg22_consensus" = "flg22")) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", face = "bold", size = 14),
        axis.line = element_line(colour = "black", 
                                 size = 0.4, linetype = "solid")) +
  scale_y_continuous(breaks = c(0,25,50,75,100),
                     limits = c(0,110)) +
  scale_fill_manual("Gram Type", values = Gram_colors) +
  geom_text(data = n_number, 
            aes(group = interaction(Gram, MAMP_Hit), y = 110, label = n), position = dodge, size = 4)


ggsave(MAMP_vs_Gram, filename = "./../Figures/Supplemental_Figure_1/Plot_MAMP_similarity_by_MAMP_and_Gram.pdf", device = cairo_pdf, width = 7, height = 3, units = "in")




#MAMP_vs_Gram_copy_number <- ggplot(subset(hold_copy_number, n != 0), aes(x = MAMP_Hit, y = n, fill = Gram)) +
#  my_ggplot_theme +
#  geom_violin(position = dodge, scale = "width") +
#  geom_boxplot(aes(group = interaction(Gram, MAMP_Hit)), position = dodge, fill = "white", color = "black", outlier.alpha = 0, width = 0.1) +
#  ylab("Number of MAMP Encoding\nProtein Sequences") +
#  scale_y_continuous(limits = c(0,14), breaks = c(0,2,4,6,8,10,12,14)) +
#  scale_x_discrete(name ="MAMP", 
#                   labels=c("csp22_consensus" = "csp22", "elf18_consensus" = "elf18",
#                            "flg22_consensus" = "flg22", 'nlp20_consensus' = "nlp20")) +
#  theme(axis.text.x = element_text(color = "black", size = 12),
#        axis.text.y = element_text(color = "black", size = 12),
#        axis.title.y = element_text(color = "black", face = "bold", size = 12),
#        axis.line = element_line(colour = "black", size = 0.4, linetype = "solid")) 

#ggsave(MAMP_vs_Gram_copy_number, filename = "./../Figures/Supplemental_Figure_1/Plot_MAMP_number_by_MAMP_and_Gram.pdf", device = cairo_pdf, width = 7, height = 2.8, units = "in")







###################################################################################################################




#-----------------------Supplemental Figure ? Plots------------------------------------------------------------------------------



##############################################
# Filter and Plot the similarity of the MAMP to the consensus seperate by bacterial class
##############################################

class_type <- list()
for (i in 1:nrow(filtered_hold_MAMP_seqs)){
  class_type[[i]] <- datasettable[datasettable$Filename %in% filtered_hold_MAMP_seqs$File_Name[i],7]
  print(class_type[[i]])
}

class_type <- unlist(class_type)

if (grepl("Class_type", colnames(filtered_hold_MAMP_seqs)) == FALSE){
  filtered_hold_MAMP_seqs <- cbind(filtered_hold_MAMP_seqs, "Class_type" = class_type)
}

dodge <- position_dodge(width = 1)
n_number <- as.data.frame(filtered_hold_MAMP_seqs %>% group_by(MAMP_Hit, Class_type) %>% summarise(n=n()))


Class_type <- ggplot(filtered_hold_MAMP_seqs, aes(x = Class_type, y = Percent_Identity, fill = MAMP_Hit)) +
  theme_classic() +
  geom_violin(position = dodge, scale = "width") +
  geom_boxplot(aes(group = interaction(Class_type, MAMP_Hit)), position = dodge, fill = "white", color = "black", outlier.alpha = 0, width = 0.1) +
  scale_fill_manual("MAMP", labels = c("csp22", "elf18", "flg22", "flgII-28", "nlp20"), values = MAMP_colors) +
  ylab("Percent AA Similarity") +
  scale_x_discrete(name ="\nBacteria Class Type", 
                   labels=c("csp22_consensus" = "csp22", 
                            "elf18_consensus" = "elf18",
                            "flg22_consensus" = "flg22",
                            'nlp20_consensus' = "nlp20")) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14, face = 'bold', family = 'Arial'),
        legend.text = element_text(color = "black", size = 10, family = 'Arial'),
        legend.title = element_text(color = "black", size = 12, face = 'bold', family = 'Arial'),
        axis.line = element_line(colour = "black", 
                                 size = 0.4, linetype = "solid")) +
  scale_y_continuous(breaks = c(0,25,50,75,100),
                     limits = c(0,110)) +
  geom_text(data = n_number, 
            aes(group = interaction(Class_type, MAMP_Hit), y = 110, label = n), position = dodge, size = 4)



ggsave(MAMP_vs_Gram_copy_number, filename = "./../Figures/Supplemental_Figure_1/Plot_MAMP_number_by_MAMP_and_Gram.pdf", device = cairo_pdf, width = 7, height = 2.8, units = "in")







# ggplot(subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus"), 
#                    aes(x = File_Name, y = Percent_Identity)) +
#  my_ggplot_theme +
#  geom_jitter(aes(color = Genera), size = 0.5) +
#geom_violin(aes(fill = MAMP_Hit), alpha = 0.9,scale = "width", trim = T) +
#geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
# ylab("Percent AA Similarity") +
#  scale_color_manual("Genera", values = Genera_colors) +
#  theme(legend.position = "none",
#        axis.text.x = element_blank(),
#        axis.text.y = element_text(color = "black", size = 14),
#        axis.title = element_text(color = "black", face = "bold", size = 14),
#        axis.line = element_line(colour = "black", 
#                                 size = 0.4, linetype = "solid")) +
#  scale_y_continuous(breaks = c(0,25,50,75,100),
#                     limits = c(0,110)) 


##############################################
# Plot the similarity of the MAMP to the consensus seperate by Gram-type, MAMP, or both - old version - ignore for now
##############################################

## NOTE! Nlp20 data is removed and is in the supplemental data 

# plots all MAMPs grouped by whether it's a gram-positive or gram-negative bacteria
n_number <- as.data.frame(filtered_hold_MAMP_seqs %>% group_by(Gram) %>% summarise(n=n()))


Gram_type <- ggplot(data = filtered_hold_MAMP_seqs, aes(x = Gram, y = Percent_Identity)) +
  my_ggplot_theme +
  geom_violin(aes(fill = Gram), alpha = 0.9, scale = "width", trim = F) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.1) +
  scale_fill_manual("Gram Type", values = Gram_colors) +
  ylab("Percent AA Similarity") +
  xlab("Gram-type") +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", size = 0.4, linetype = "solid")) +
  scale_y_continuous(breaks = c(0,25,50,75,100),
                     limits = c(0,110)) +
  geom_text(data = n_number, aes(x = Gram, y = 110, label = n), size = 4)



Gram_type <- Gram_type + 
  ggplot(filtered_hold_MAMP_seqs, aes(color = Gram, x = Percent_Identity)) + 
  stat_ecdf(geom = "step", size = 1.25) +
  my_ggplot_theme +
  ylab("Cumulative Probability") +
  xlab("Percent AA Similarity") +
  scale_color_manual("Gram Type", values = Gram_colors) +
  theme(legend.position = "none")


# save plot
ggsave(Gram_type, filename = "./../Figures/Supplemental_Figure_1/Plot_MAMP_similarity_by_Gram_type.pdf", device = cairo_pdf, width = 5.5, height = 2.8, units = "in")





############################################################################ 
# repeat on subset of data -> MAMP has had to occur at teast 10 times -- IGNORE
############################################################################ 


parising_MAMP_all_by_all_comparisons <- function(MAMP_df_in, occurance_number){
  
}

csp22_above_10 <- subset(csp22_occurance, csp22_occurance$n > 10)

sim_to_csp22_consensus <- list()
for (i in 1:nrow(csp22_above_10)){
  test2 <- Biostrings::pairwiseAlignment(consensus_csp22, csp22_above_10$MAMP_Sequence[[i]], type = "global")
  sim_to_csp22_consensus[[i]] <- Biostrings::pid(test2, type = "PID1")
}

csp22_above_10 <- cbind(csp22_above_10, "Similarity_to_con_cso22" = unlist(sim_to_csp22_consensus))


csp22_above_10 <- csp22_above_10[sort(csp22_above_10$Similarity_to_con_cso22, decreasing = T, index.return = T)$ix,]


# initiate matrix
empty_csp22_matrix_above10 <- matrix(0, 
                                     nrow = length(csp22_above_10$MAMP_Sequence), 
                                     ncol = length(csp22_above_10$MAMP_Sequence),
                                     dimnames = list(csp22_above_10$MAMP_Sequence, csp22_above_10$MAMP_Sequence))

subset_sim_score <- csp22_sim_score[csp22_sim_score$MAMP_v1 %in% csp22_above_10$MAMP_Sequence,]
subset_sim_score <- subset_sim_score[subset_sim_score$MAMP_v2 %in% csp22_above_10$MAMP_Sequence,]

# fill in upper triangle
empty_csp22_matrix_above10[as.matrix(subset_sim_score[c(1,2)])] <- subset_sim_score$Percent_sim


# subset eptitope versions with more than 10 occruanecs
csp22_above_10_by_genera <- rbind(csp22_occurance_by_genera %>% filter(n.Streptomyces > 10), 
                                  csp22_occurance_by_genera %>% filter(n.Rhodococcus > 10),
                                  csp22_occurance_by_genera %>% filter(n.Clavibacter > 10),
                                  csp22_occurance_by_genera %>% filter(n.Xanthomonas > 10),
                                  csp22_occurance_by_genera %>% filter(n.Leifsonia > 10),
                                  csp22_occurance_by_genera %>% filter(n.Agrobacterium > 10),
                                  csp22_occurance_by_genera %>% filter(n.Curtobacterium > 10),
                                  csp22_occurance_by_genera %>% filter(n.Ralstonia > 10),
                                  csp22_occurance_by_genera %>% filter(n.Rathayibacter > 10),
                                  csp22_occurance_by_genera %>% filter(n.Pseudomonas > 10),
                                  csp22_occurance_by_genera %>% filter(n.Pectobacterium > 10),
                                  csp22_occurance_by_genera %>% filter(n.Dickeya > 10),
                                  csp22_occurance_by_genera %>% filter(n.Erwinia >10)
)


# reorder subset to match same order as MAMP occurance
csp22_above_10_by_genera <- csp22_above_10_by_genera[match(csp22_above_10$MAMP_Sequence, csp22_above_10_by_genera$MAMP_Sequence),]


ha <- rowAnnotation("MAMP \noccurance" = anno_barplot(csp22_above_10$n, width = unit(2, "cm")),
                    Similarity_to_con = (csp22_above_10$Similarity_to_con_cso22),
                    "Streptomyces" = csp22_above_10_by_genera$n.Streptomyces,
                    "Rhodococcus" = csp22_above_10_by_genera$n.Rhodococcus,
                    "Clavibacter" = csp22_above_10_by_genera$n.Clavibacter,
                    "Leifsonia" = csp22_above_10_by_genera$n.Leifsonia,
                    "Curtobacterium" = csp22_above_10_by_genera$n.Curtobacterium,
                    "Rathayibacter" = csp22_above_10_by_genera$n.Rathayibacter,
                    
                    
                    "Xanthomonas" = csp22_above_10_by_genera$n.Xanthomonas,
                    "Agrobacterium" = csp22_above_10_by_genera$n.Agrobacterium,
                    "Ralstonia" = csp22_above_10_by_genera$n.Ralstonia,
                    "Pseudomonas" = csp22_above_10_by_genera$n.Pseudomonas,
                    
                    "Pectobacterium" = csp22_above_10_by_genera$n.Pectobacterium,
                    "Dickeya" = csp22_above_10_by_genera$n.Dickeya,
                    "Erwinia" = csp22_above_10_by_genera$n.Erwinia,
                    
                    
                    gp = gpar(col = "black", size = 1)
                    
)

ComplexHeatmap::Heatmap(empty_csp22_matrix_above10,
                        #row_km = 5,
                        #row_km_repeats = 1000,
                        #column_km = ,
                        
                        #column_split = 4,
                        right_annotation = ha,
                        border = T,
                        
                        show_column_dend = F,
                        row_dend_width = unit(4.5,"cm"),                       
                        show_row_names = T,
                        show_column_names = F, 
                        row_names_gp = gpar(fontsize = 7),
                        
                        row_dend_reorder = T,
                        column_dend_reorder = T)




############################################################################ 
# repeat on subset of data -> MAMP has had to occur at teast 30 times
############################################################################ 




csp22_above_30 <- subset(csp22_occurance, csp22_occurance$n > 50)

sim_to_csp22_consensus <- list()
for (i in 1:nrow(csp22_above_30)){
  test2 <- Biostrings::pairwiseAlignment(consensus_csp22, csp22_above_30$MAMP_Sequence[[i]], type = "global")
  sim_to_csp22_consensus[[i]] <- Biostrings::pid(test2, type = "PID1")
}

csp22_above_30 <- cbind(csp22_above_30, "Similarity_to_con_cso22" = unlist(sim_to_csp22_consensus))


csp22_above_30 <- csp22_above_30[sort(csp22_above_30$Similarity_to_con_cso22, decreasing = T, index.return = T)$ix,]


# initiate matrix
empty_csp22_matrix_above30 <- matrix(0, 
                                     nrow = length(csp22_above_30$MAMP_Sequence), 
                                     ncol = length(csp22_above_30$MAMP_Sequence),
                                     dimnames = list(csp22_above_30$MAMP_Sequence, csp22_above_30$MAMP_Sequence))

subset_sim_score <- csp22_sim_score[csp22_sim_score$MAMP_v1 %in% csp22_above_30$MAMP_Sequence,]
subset_sim_score <- subset_sim_score[subset_sim_score$MAMP_v2 %in% csp22_above_30$MAMP_Sequence,]

# fill in upper triangle
empty_csp22_matrix_above30[as.matrix(subset_sim_score[c(1,2)])] <- subset_sim_score$Percent_sim


# subset eptitope versions with more than 10 occruanecs
csp22_above_30_by_genera <- rbind(csp22_occurance_by_genera %>% filter(n.Streptomyces > 10), 
                                  csp22_occurance_by_genera %>% filter(n.Rhodococcus > 10),
                                  csp22_occurance_by_genera %>% filter(n.Clavibacter > 10),
                                  csp22_occurance_by_genera %>% filter(n.Xanthomonas > 10),
                                  csp22_occurance_by_genera %>% filter(n.Leifsonia > 10),
                                  csp22_occurance_by_genera %>% filter(n.Agrobacterium > 10),
                                  csp22_occurance_by_genera %>% filter(n.Curtobacterium > 10),
                                  csp22_occurance_by_genera %>% filter(n.Ralstonia > 10),
                                  csp22_occurance_by_genera %>% filter(n.Rathayibacter > 10),
                                  csp22_occurance_by_genera %>% filter(n.Pseudomonas > 10)
)


# reorder subset to match same order as MAMP occurance
csp22_above_30_by_genera <- csp22_above_30_by_genera[match(csp22_above_30$MAMP_Sequence, csp22_above_30_by_genera$MAMP_Sequence),]


ha <- rowAnnotation("MAMP \noccurance" = anno_barplot(csp22_above_30$n, width = unit(2, "cm")),
                    Similarity_to_con = (csp22_above_30$Similarity_to_con_cso22),
                    "Streptomyces" = csp22_above_30_by_genera$n.Streptomyces,
                    "Rhodococcus" = (csp22_above_30_by_genera$n.Rhodococcus),
                    "Clavibacter" = (csp22_above_30_by_genera$n.Clavibacter),
                    "Xanthomonas" = (csp22_above_30_by_genera$n.Xanthomonas),
                    "Leifsonia" = (csp22_above_30_by_genera$n.Leifsonia),
                    "Agrobacterium" = (csp22_above_30_by_genera$n.Agrobacterium),
                    "Curtobacterium" = (csp22_above_30_by_genera$n.Curtobacterium),
                    "Ralstonia" = (csp22_above_30_by_genera$n.Ralstonia),
                    "Rathayibacter" = (csp22_above_30_by_genera$n.Rathayibacter),
                    "Pseudomonas" = (csp22_above_30_by_genera$n.Pseudomonas)
)


color_species <- c("Leifsonia" = "#ff94af",
                   "Clavibacter" = "#ffd65c",
                   "Curtobacterium" = "#b589d6",
                   "Rhodococcus" = "#73bfe6",
                   "Rathayibacter" = "#9e9e9e",
                   "Streptomyces" = "#7bc98f",
                   "Ralstonia" = "#ab234c",
                   "Xanthomonas" =  "#e3534c",
                   "Pseudomonas" = "#026178",
                   "Agrobacterium" = "#507B00"
)

ComplexHeatmap::Heatmap(empty_csp22_matrix_above30,
                        #row_km = 5,
                        #row_km_repeats = 1000,
                        #column_km = ,
                        
                        #column_split = 4,
                        right_annotation = ha,
                        
                        
                        show_column_dend = F,
                        row_dend_width = unit(4.5,"cm"),                       
                        show_row_names = T,
                        show_column_names = F, 
                        row_names_gp = gpar(fontsize = 10),
                        
                        row_dend_reorder = T,
                        column_dend_reorder = T)





############################################################################

test_weblogo <- c("ANGTVKWFNAEKGFGFITVDGG", "ATGTVKWFNAEKGFGFIAQDGG", "ATGTVKWFNAEKGFGFIEQDGG", "ATGTVKWFNSEKGFGFIEQDGG",
                  "ASGTVKWFNAEKGFGFIEQDGG", "AAGTVKWFNAEKGFGFIEQDGG", "ATGTVKWFNAEKGFGFIAQEGG", "AQGTVKWFNAEKGFGFIAPEDG", 
                  "ASGTVKWFNSEKGFGFIEQEGG", "TQGSVKWFNGEKGFGFIEQDGG", "ANGTVKWFNGEKGFGFITVDAV", "ANGTVKWFNAEKGYGFITVDGG",
                  "AQGTVKWFNAEKGYGFIAVDGG", "ANGTVKWFNAEKGYGFITVDGS")

test_weblogo2 <- c("PTGKVKFYDDEKGFGFISTDDG", "PTGKVKFYDDEKGFGFISSDDG", "PTGKVKFYDDQKGFGFISGDDG", "PTGKVKFYDDDKGFGFITGDDG",
                   "PTGKVKFYDDQKGFGFITGDDG", "PTGKVKFYDEEKGFGFISTDDG", "PTGKVKFYDEEKGFGFISSDDG", "PTGKVKFYDEEKGFGFISTDEG", 
                   "PTGKVKFYDEDKGFGFISSDDG", "PTGKVKWYDVDKGFGFLSQEEG")



##############################################
# how often does each version of the epitope occur
##############################################

elf18_occurance <- as.data.frame(elf18_MAMP_seqs %>% group_by(MAMP_Sequence) %>% summarise(n=n()))
csp22_occurance <- as.data.frame(csp22_MAMP_seqs %>% group_by(MAMP_Sequence) %>% summarise(n=n()))
flg22_occurance <- as.data.frame(flg22_MAMP_seqs %>% group_by(MAMP_Sequence) %>% summarise(n=n()))








# Comparing 


##############################################
# consensus MAMP
##############################################

#consensus_csp22 <- c("AVGTVKWFNAEKGFGFITPDDG")
#consensus_flg22 <- c("QRLSTGSRINSAKDDAAGLQIA")
#consensus_elf18 <- c("SKEKFERTKPHVNVGTIG")



##############################################
# make weblogos of all variants based - IGNORE
##############################################


### plotting elf18 variants as weblogo to show changes in a position basis - part of Figure 1G
make_me_a_weblogo(unique(subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "elf18_consensus")[[4]])) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))


### plotting csp22 variants as weblogo to show changes in a position basis - part of Figure 1G
make_me_a_weblogo(unique(subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")[[4]])) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))


### plotting variants as weblogo to show changes in a position basis  
make_me_a_weblogo(unique(subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flg22_consensus")[[4]])) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10)) 


make_me_a_weblogo(unique(subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flgII-28")[[4]])) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))



##############################################
# make all-by-all comparisons of flg22 variants - IGNORE THIS FOR NOW
##############################################



### summary of elf18 heatmap as a violin plot - part of Figure 1G
### remove comparisons to itself (i.e artifical 100% score)
ggplot(elf18_sim_score %>% filter(MAMP_v1!=MAMP_v2), aes(x = "", y = Percent_sim)) + 
  geom_violin(fill = "grey", scale = "width", trim = T) +
  my_ggplot_theme +
  ylab("Percent Similarity") +
  xlab("") +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10)) 


### summary of csp22 heatmap as a violin plot - part of Figure 1G
### remove comparisons to itself (i.e artifical 100% score)
ggplot(csp22_sim_score %>% filter(MAMP_v1!=MAMP_v2), aes(x = "", y = Percent_sim)) + 
  geom_violin(fill = "grey", scale = "width", trim = T) +
  my_ggplot_theme +
  ylab("Percent Similarity") +
  xlab("") +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))


ggplot(flg22_sim_score %>% filter(MAMP_v1!=MAMP_v2), aes(x = "", y = Percent_sim)) + 
  geom_violin(fill = "grey", scale = "width", trim = T) +
  my_ggplot_theme +
  ylab("Percent Similarity") +
  xlab("") +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))



ggplot(flgII_28_sim_score %>% filter(MAMP_v1!=MAMP_v2), aes(x = "", y = Percent_sim)) + 
  geom_violin(fill = "grey", scale = "width", trim = T) +
  my_ggplot_theme +
  ylab("Percent Similarity") +
  xlab("") +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  theme(axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))





##############################################
# Load colors - Set tip labels to match other figures
##############################################


######### csp22 peptide
clavibacter_csps <- subset(All_target_by_annotation, All_target_by_annotation$MAMP_Hit == "csp22_consensus")
clavibacter_csps <- subset(clavibacter_csps, clavibacter_csps$Genera == "Clavibacter")


# 
pb <- txtProgressBar(min = 0, max = nrow(clavibacter_csps), style = 3)
clavibacter_sim_score <- data.frame("Protein_1" = character(0), "Protein_2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:nrow(clavibacter_csps)){
  for (j in 2:nrow(clavibacter_csps)){
    test2 <- Biostrings::pairwiseAlignment(clavibacter_csps$seq[[i]], clavibacter_csps$seq[[j]], type = "global", substitutionMatrix = BLOSUM62)
    hold_score <- Biostrings::pid(test2, type = "PID1")
    clavibacter_sim_score <- rbind(clavibacter_sim_score, data.frame("Protein_1" = paste(clavibacter_csps$Filename[[i]], clavibacter_csps$Protein_Name[[i]], sep = "|"),
                                                                     "Protein_2" = paste(clavibacter_csps$Filename[[j]], clavibacter_csps$Protein_Name[[j]], sep = "|"),
                                                                     "Percent_sim" = hold_score))
  }
  setTxtProgressBar(pb, i)
} 



clavibacter_csp_names <- sort(unique(as.character(unlist(clavibacter_sim_score[1:2]))))

empty_clavibacter_csp22_matrix <- matrix(0, 
                                         nrow = length(clavibacter_csp_names), 
                                         ncol = length(clavibacter_csp_names),
                                         dimnames = list(clavibacter_csp_names, clavibacter_csp_names))

# fill in upper triangle
empty_clavibacter_csp22_matrix[as.matrix(clavibacter_sim_score[c(1,2)])] <- clavibacter_sim_score$Percent_sim



clav_test <- ComplexHeatmap::Heatmap(empty_clavibacter_csp22_matrix,
                                     #row_km = 5,
                                     #row_km_repeats = 2000,
                                     #column_km = ,
                                     
                                     #column_split = 4,
                                     
                                     show_column_dend = F,
                                     row_dend_width = unit(4.5,"cm"),                       
                                     show_row_names = T,
                                     show_column_names = F, 
                                     row_names_gp = gpar(fontsize = 2.5),
                                     border = T,
                                     row_dend_reorder = T,
                                     column_dend_reorder = T)



hold_order <- row_order(clav_test)


All_protein_values <- data.frame("Protein_1" = character(0), "Protein_2" = character(0), 
                                 "Percent_Similarity" = numeric(0), "Cluster" = character(0))

for (i in 1:length(hold_order)){
  
  cluster_1 <- empty_clavibacter_csp22_matrix[hold_order[[i]],hold_order[[i]]]
  cluster_1_df <- do.call(rbind, strsplit(colnames(cluster_1), split = "|", fixed = T))
  cluster_1 <- as.data.frame(matrix(cluster_1, dimnames=list(t(outer(colnames(cluster_1), rownames(cluster_1), FUN=paste)), NULL)))
  cluster_1 <- cbind(cluster_1, "Cluster" = rep("Cluster_1", nrow(cluster_1)))
  
  cluster_2 <- empty_clavibacter_csp22_matrix[hold_order$`2`,hold_order$`2`]
  cluster_2_df <- do.call(rbind, strsplit(colnames(cluster_2), split = "|", fixed = T))
  cluster_2 <- as.data.frame(matrix(cluster_2, dimnames=list(t(outer(colnames(cluster_2), rownames(cluster_2), FUN=paste)), NULL)))
  cluster_2 <- cbind(cluster_2, "Cluster" = rep("Cluster_2", nrow(cluster_2)))
  
  cluster_3 <- empty_clavibacter_csp22_matrix[hold_order$`3`,hold_order$`3`]
  cluster_3_df <- do.call(rbind, strsplit(colnames(cluster_3), split = "|", fixed = T))
  cluster_3 <- as.data.frame(matrix(cluster_3, dimnames=list(t(outer(colnames(cluster_3), rownames(cluster_3), FUN=paste)), NULL)))
  cluster_3 <- cbind(cluster_3, "Cluster" = rep("Cluster_3", nrow(cluster_3)))
  
  All_protein_values <- rbind(cluster_1, cluster_2, cluster_3) 
  ggplot(All_protein_values, aes(x=Cluster, y=V1)) + geom_boxplot() +
    my_ggplot_theme
  
  colnames(cluster_1_df) <- c("Filename", "Protein_Name")
  cluster_1_df <- as.data.frame(cluster_1_df)
  
  
  
  cluster_df_1 <- clavibacter_csps[clavibacter_csps$Filename %in% cluster_1_df$Filename &
                                     clavibacter_csps$Protein_Name %in% cluster_1_df$Protein_Name,]
  
  
  
  ##############################################
  # Load colors - Set tip labels to match other figures
  ##############################################
  
  
  ######### csp22 peptide
  agrobacterium_csps <- subset(All_target_by_annotation, All_target_by_annotation$MAMP_Hit == "csp22_consensus")
  agrobacterium_csps <- subset(agrobacterium_csps, agrobacterium_csps$Genera == "Rhodococcus")
  
  
  # 
  pb <- txtProgressBar(min = 0, max = nrow(agrobacterium_csps), style = 3)
  agrobacterium_sim_score <- data.frame("Protein_1" = character(0), "Protein_2" = character(0), "Percent_sim" = numeric(0))
  for (i in 1:nrow(agrobacterium_csps)){
    for (j in 2:nrow(agrobacterium_csps)){
      test2 <- Biostrings::pairwiseAlignment(agrobacterium_csps$seq[[i]], agrobacterium_csps$seq[[j]], type = "global", substitutionMatrix = BLOSUM62)
      hold_score <- Biostrings::pid(test2, type = "PID1")
      agrobacterium_sim_score <- rbind(agrobacterium_sim_score, data.frame("Protein_1" = paste(agrobacterium_csps$Filename[[i]], agrobacterium_csps$Protein_Name[[i]], sep = "|"),
                                                                           "Protein_2" = paste(agrobacterium_csps$Filename[[j]], agrobacterium_csps$Protein_Name[[j]], sep = "|"),
                                                                           "Percent_sim" = hold_score))
    }
    setTxtProgressBar(pb, i)
  } 
  
  
  
  agrobacterium_csp_names <- sort(unique(as.character(unlist(agrobacterium_sim_score[1:2]))))
  
  empty_agrobacterium_csp22_matrix <- matrix(0, 
                                             nrow = length(agrobacterium_csp_names), 
                                             ncol = length(agrobacterium_csp_names),
                                             dimnames = list(agrobacterium_csp_names, agrobacterium_csp_names))
  
  # fill in upper triangle
  empty_agrobacterium_csp22_matrix[as.matrix(agrobacterium_sim_score[c(1,2)])] <- agrobacterium_sim_score$Percent_sim
  
  
  
  ComplexHeatmap::Heatmap(empty_agrobacterium_csp22_matrix,
                          #row_km = 5,
                          #row_km_repeats = 1000,
                          #column_km = ,
                          
                          #column_split = 4,
                          
                          show_column_dend = F,
                          row_dend_width = unit(4.5,"cm"),                       
                          show_row_names = T,
                          show_column_names = F, 
                          row_names_gp = gpar(fontsize = 2.5),
                          border = T,
                          row_dend_reorder = T,
                          column_dend_reorder = T)
  
