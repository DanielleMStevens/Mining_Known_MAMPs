#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 11/12/20222
# Script Purpose: Figure 1G
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------

#-----------------------Figure 2 Plots------------------------------------------------------------------------------

##############################################
# create a list of each MAMP and determine the unique variants for comparison
##############################################


elf18_MAMP_seqs <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "elf18_consensus")
elf18_MAMP_seqs <- subset(elf18_MAMP_seqs, elf18_MAMP_seqs$Percent_Identity > 50)
elf18_unique <- unique(elf18_MAMP_seqs$MAMP_Sequence)

csp22_MAMP_seqs <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")
csp22_MAMP_seqs <- subset(csp22_MAMP_seqs, csp22_MAMP_seqs$Percent_Identity > 20)
csp22_unique <- unique(csp22_MAMP_seqs$MAMP_Sequence)

flg22_MAMP_seqs <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flg22_consensus")
flg22_unique <- unique(flg22_occurance$MAMP_Sequence)

flgII_28_MAMP_seqs <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "flgII-28")
flgII_28_unique <- unique(flgII_28_MAMP_seqs$MAMP_Sequence)

nlp20_MAMP_seqs <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "nlp20_consensus")
nlp20_unique <- unique(nlp20_MAMP_seqs$MAMP_Sequence)

##############################################
# make all-by-all comparisons of each MAMP and their varaints
##############################################


### Compare all by all local similarity score - elf18 variant comparisons
pb <- txtProgressBar(min = 0, max = length(elf18_unique), style = 3)
elf18_sim_score <- data.frame("MAMP_v1" = character(0), "MAMP_v2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:length(elf18_unique)){
  for (j in 1:length(elf18_unique)){
    elf18_sim_score <- rbind(elf18_sim_score, data.frame("MAMP_v1" = elf18_unique[[i]],
                                                         "MAMP_v2" = elf18_unique[[j]],
                                                         "Percent_sim" = Biostrings::pid(Biostrings::pairwiseAlignment(elf18_unique[[i]], elf18_unique[[j]], 
                                                                                                                       type = "global", substitutionMatrix = BLOSUM62), type = "PID1")))
  }
  setTxtProgressBar(pb, i)
  #print(i) - for debugging 
} 


### Compare all by all local similarity score - csp22 variant comparisons
### this takes a long time to complete -1.5-2 hours
pb <- txtProgressBar(min = 0, max = length(csp22_unique), style = 3)
csp22_sim_score <- data.frame("MAMP_v1" = character(0), "MAMP_v2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:length(csp22_unique)){
  for (j in 1:length(csp22_unique)){
    csp22_sim_score <- rbind(csp22_sim_score, data.frame("MAMP_v1" = csp22_unique[[i]],
                                                         "MAMP_v2" = csp22_unique[[j]],
                                                         "Percent_sim" = Biostrings::pid(Biostrings::pairwiseAlignment(csp22_unique[[i]], csp22_unique[[j]], 
                                                                                                                       type = "global", substitutionMatrix = BLOSUM62), type = "PID1")))
  }
  setTxtProgressBar(pb, i)
} 


### Compare all by all local similarity score - flg22 variant comparisons
pb <- txtProgressBar(min = 0, max = length(flg22_unique), style = 3)
flg22_sim_score <- data.frame("MAMP_v1" = character(0), "MAMP_v2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:length(flg22_unique)){
  for (j in 1:length(flg22_unique)){
    flg22_sim_score <- rbind(flg22_sim_score, data.frame("MAMP_v1" = flg22_unique[[i]],
                                                         "MAMP_v2" = flg22_unique[[j]],
                                                         "Percent_sim" =  Biostrings::pid(Biostrings::pairwiseAlignment(flg22_unique[[i]], flg22_unique[[j]], 
                                                                                                                        type = "global", substitutionMatrix = BLOSUM62), type = "PID1") ))
  }
  setTxtProgressBar(pb, i)
} 


### Compare all by all local similarity score - flgII-28 variant comparisons
pb <- txtProgressBar(min = 0, max = length(flgII_28_unique), style = 3)
flgII_28_sim_score <- data.frame("MAMP_v1" = character(0), "MAMP_v2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:length(flgII_28_unique)){
  for (j in 1:length(flgII_28_unique)){
    flgII_28_sim_score <- rbind(flgII_28_sim_score, data.frame("MAMP_v1" = flgII_28_unique[[i]],
                                                               "MAMP_v2" = flgII_28_unique[[j]],
                                                               "Percent_sim" =  Biostrings::pid(Biostrings::pairwiseAlignment(flgII_28_unique[[i]], flgII_28_unique[[j]], 
                                                                                                                        type = "global", substitutionMatrix = BLOSUM62), type = "PID1") ))
  }
  setTxtProgressBar(pb, i)
} 

### Compare all by all local similarity score - nlp20 variant comparisons
pb <- txtProgressBar(min = 0, max = length(nlp20_unique), style = 3)
nlp20_sim_score <- data.frame("MAMP_v1" = character(0), "MAMP_v2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:length(nlp20_unique)){
  for (j in 1:length(nlp20_unique)){
    nlp20_sim_score <- rbind(nlp20_sim_score, data.frame("MAMP_v1" = nlp20_unique[[i]],
                                                         "MAMP_v2" = nlp20_unique[[j]],
                                                         "Percent_sim" = Biostrings::pid(Biostrings::pairwiseAlignment(nlp20_unique[[i]], nlp20_unique[[j]], 
                                                                                                                       type = "global", substitutionMatrix = BLOSUM62), type = "PID1")))
  }
  setTxtProgressBar(pb, i)
  #print(i) - for debugging 
} 

##############################################
### reorganize table for plotting 
##############################################


### reorganize table for plotting - elf18
elf18_Names <- sort(unique(as.character(unlist(elf18_sim_score[1]))))
empty_elf18_matrix <- matrix(0, nrow = length(elf18_Names), ncol = length(elf18_Names),
                             dimnames = list(elf18_Names, elf18_Names))


### reorganize table for plotting - csp22
csp22_Names <- sort(unique(as.character(unlist(csp22_sim_score[1]))))
empty_csp22_matrix <- matrix(0, nrow = length(csp22_Names), ncol = length(csp22_Names),
                             dimnames = list(csp22_Names, csp22_Names))


### reorganize table for plotting - flg22
flg22_Names <- sort(unique(as.character(unlist(flg22_sim_score[1]))))
empty_flg22_matrix <- matrix(0, nrow = length(flg22_Names), ncol = length(flg22_Names),
                             dimnames = list(flg22_Names, flg22_Names))


### reorganize table for plotting - flgII=28
flgII_28_Names <- sort(unique(as.character(unlist(flgII_28_sim_score[1]))))
empty_flgII_28_matrix <- matrix(0, nrow = length(flgII_28_Names), ncol = length(flgII_28_Names),
                             dimnames = list(flgII_28_Names, flgII_28_Names))

### reorganize table for plotting - nlp20
nlp20_Names <- sort(unique(as.character(unlist(nlp20_sim_score[1]))))
empty_nlp20_matrix <- matrix(0, nrow = length(nlp20_Names), ncol = length(nlp20_Names),
                                dimnames = list(nlp20_Names, nlp20_Names))

### fill in upper triangle
empty_elf18_matrix[as.matrix(elf18_sim_score[1:2])] <- elf18_sim_score$Percent_sim
empty_csp22_matrix[as.matrix(csp22_sim_score[1:2])] <- csp22_sim_score$Percent_sim
empty_flg22_matrix[as.matrix(flg22_sim_score[1:2])] <- flg22_sim_score$Percent_sim
empty_flgII_28_matrix[as.matrix(flgII_28_sim_score[1:2])] <- flgII_28_sim_score$Percent_sim
empty_nlp20_matrix[as.matrix(nlp20_sim_score[1:2])] <- nlp20_sim_score$Percent_sim



##############################################
# plotting all-by-all comparisons
##############################################

### elf18 variants 
col_fun = colorRamp2(c( 0, 50, 100), c("#00008b", "white", "#8b0000"))


elf18_heatmap <- ComplexHeatmap::Heatmap(empty_elf18_matrix,
                                         col = col_fun,
                                         show_column_dend = F,
                                         show_row_dend = F,
                                         show_row_names = F,
                                         show_column_names = F, 
                                         border = T,
                                         row_dend_reorder = T,
                                         column_dend_reorder = T)


csp22_heatmap <- ComplexHeatmap::Heatmap(empty_csp22_matrix,
                                         col = col_fun,
                                         show_column_dend = F,
                                         show_row_dend = F,
                                         show_row_names = F,
                                         show_column_names = F, 
                                         border = T,
                                         row_dend_reorder = T,
                                         column_dend_reorder = T)


flg22_heatmap <- ComplexHeatmap::Heatmap(empty_flg22_matrix,
                                         col = col_fun,
                                         show_column_dend = F,
                                         show_row_dend = F,
                                         show_row_names = F,
                                         show_column_names = F, 
                                         border = T,
                                         row_dend_reorder = T,
                                         column_dend_reorder = T)


flgII_28_heatmap <- ComplexHeatmap::Heatmap(empty_flgII_28_matrix,
                                            col = col_fun,
                                            show_column_dend = F,
                                            show_row_dend = F,
                                            show_row_names = F,
                                            show_column_names = F, 
                                            border = T,
                                            row_dend_reorder = T,
                                            column_dend_reorder = T)

nlp20_heatmap <- ComplexHeatmap::Heatmap(empty_nlp20_matrix,
                                            col = col_fun,
                                            show_column_dend = F,
                                            show_row_dend = F,
                                            show_row_names = F,
                                            show_column_names = F, 
                                            border = T,
                                            row_dend_reorder = T,
                                            column_dend_reorder = T)



### ------------- elf18 heatmap as part of Figure 1G -- export as 2.5 x 2.2 inch pdf     
draw(elf18_heatmap, show_heatmap_legend = FALSE, padding = unit(c(7, 7, 7, 7), "mm"))
draw(csp22_heatmap, show_heatmap_legend = FALSE, padding = unit(c(7, 7, 7, 7), "mm"))
draw(flg22_heatmap, show_heatmap_legend = FALSE, padding = unit(c(5, 5, 5, 5), "mm"))
draw(flgII_28_heatmap, show_heatmap_legend = FALSE, padding = unit(c(5, 5, 5, 5), "mm"))
draw(nlp20_heatmap, show_heatmap_legend = FALSE, padding = unit(c(5, 5, 5, 5), "mm"))




##############################################
# View all-by-all as a density plot
##############################################


# Density plots
(ggplot(csp22_sim_score, aes(x=Percent_sim)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=mean(Percent_sim)), color="blue",
             linetype="dashed")+
  labs(title="csp22",x="Percent Similarity", y = "Density")+
  theme_classic() +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,20,40,60,80,100)) ) /


(ggplot(elf18_sim_score, aes(x=Percent_sim)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=mean(Percent_sim)), color="blue",
             linetype="dashed")+
  labs(title="elf18",x="Percent Similarity", y = "Density")+
  theme_classic() +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,20,40,60,80,100)) ) /


  (ggplot(flg22_sim_score, aes(x=Percent_sim)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=mean(Percent_sim)), color="blue",
             linetype="dashed")+
  labs(title="flg22",x="Percent Similarity", y = "Density")+
  theme_classic() +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,20,40,60,80,100)) ) /


(ggplot(flgII_28_sim_score, aes(x=Percent_sim)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=mean(Percent_sim)), color="blue",
             linetype="dashed")+
  labs(title="flgII-28",x="Percent Similarity", y = "Density")+
  theme_classic() +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,20,40,60,80,100))) / 
 
  ggplot(nlp20_sim_score, aes(x=Percent_sim)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=mean(Percent_sim)), color="blue",
             linetype="dashed")+
  labs(title="nlp20",x="Percent Similarity", y = "Density")+
  theme_classic() +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0,20,40,60,80,100))
 


  

