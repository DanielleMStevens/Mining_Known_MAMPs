#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 11/12/20222
# Script Purpose: Figure 1G
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------

#-----------------------Figure 1 Plots------------------------------------------------------------------------------

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

### fill in upper triangle
empty_elf18_matrix[as.matrix(elf18_sim_score[1:2])] <- elf18_sim_score$Percent_sim
empty_csp22_matrix[as.matrix(csp22_sim_score[1:2])] <- csp22_sim_score$Percent_sim
empty_flg22_matrix[as.matrix(flg22_sim_score[1:2])] <- flg22_sim_score$Percent_sim
empty_flgII_28_matrix[as.matrix(flgII_28_sim_score[1:2])] <- flgII_28_sim_score$Percent_sim



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



### ------------- elf18 heatmap as part of Figure 1G -- export as 2.5 x 2.2 inch pdf     
draw(elf18_heatmap, show_heatmap_legend = FALSE, padding = unit(c(7, 7, 7, 7), "mm"))
draw(csp22_heatmap, show_heatmap_legend = FALSE, padding = unit(c(7, 7, 7, 7), "mm"))
draw(flg22_heatmap, show_heatmap_legend = FALSE, padding = unit(c(5, 5, 5, 5), "mm"))
draw(flgII_28_heatmap, show_heatmap_legend = FALSE, padding = unit(c(5, 5, 5, 5), "mm"))



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


#########################################################################
# reorganize table for plotting - csp22
#########################################################################



#ha <- rowAnnotation("MAMP \nOccurance" = anno_barplot(elf18_occurance$n, width = unit(1, "cm")))




#####################
(ggplot(elf18_sim_score, aes(x = "", y = Percent_sim)) + geom_violin(fill = "grey", alpha = 0.7, trim = T) + 
   scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105)) + my_ggplot_theme) +
(ggplot(csp22_sim_score, aes(x = "", y = Percent_sim)) + geom_violin(fill = "grey", alpha = 0.7, trim = T) + 
   scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105)) + my_ggplot_theme) +
(ggplot(flg22_sim_score, aes(x = "", y = Percent_sim)) + geom_violin(fill = "grey", alpha = 0.7, trim = T) + 
   scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105)) + my_ggplot_theme) 
  
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




###################################################################################################

#hold_cspD <- c("MQNGKVKWFNNEKGFGFIEVEGGDDVFVHFTAIEGDGYKSLEEGQEVSFEIVEGNRGPQASNVVKL")
#hold_cspB <- c("MLEGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTKEA")

test2 <- Biostrings::pairwiseAlignment(hold_cspD, All_target_by_annotation$seq[4], type = "global")
test2
Biostrings::pid(test2, type = "PID1")








# Comparing 


##############################################
# consensus MAMP
##############################################

#consensus_csp22 <- c("AVGTVKWFNAEKGFGFITPDDG")
#consensus_flg22 <- c("QRLSTGSRINSAKDDAAGLQIA")
#consensus_elf18 <- c("SKEKFERTKPHVNVGTIG")

