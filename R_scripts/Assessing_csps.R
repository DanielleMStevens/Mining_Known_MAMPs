#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------




##############################################
# Load colors - Set tip labels to match other figures
##############################################

csp22_MAMP_seqs <- subset(hold_MAMP_seqs, hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")

# which csp22 sequences are unique
csp22_unique <- unique(csp22_MAMP_seqs$MAMP_Sequence)


# this takes a long time to complete -1.5-2 hours
pb <- txtProgressBar(min = 0, max = length(csp22_unique), style = 3)
csp22_sim_score <- data.frame("MAMP_v1" = character(0), "MAMP_v2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:length(csp22_unique)){
  for (j in 2:length(csp22_unique)){
    test2 <- Biostrings::pairwiseAlignment(csp22_unique[[i]], csp22_unique[[j]], type = "global")
    hold_score <- Biostrings::pid(test2, type = "PID1")
    csp22_sim_score <- rbind(csp22_sim_score, data.frame("MAMP_v1" = csp22_unique[[i]],
                                                         "MAMP_v2" = csp22_unique[[j]],
                                                         "Percent_sim" = hold_score))
  }
  setTxtProgressBar(pb, i)
  print(i)
} 


##############################################
# can I order based on simialrity from consnensus
##############################################

consensus_csp22 <- c("AVGTVKWFNAEKGFGFITPDDG")



#########################################################################
# reorganize table for plotting
#########################################################################

csp22_Names <- sort(unique(as.character(unlist(csp22_sim_score[1:2]))))

empty_csp22_matrix <- matrix(0, 
                             nrow = length(csp22_Names), 
                             ncol = length(csp22_Names),
                             dimnames = list(csp22_Names, csp22_Names))

# fill in upper triangle
empty_csp22_matrix[as.matrix(csp22_sim_score[c(1,2)])] <- csp22_sim_score$Percent_sim

csp22_occurance <- as.data.frame(csp22_MAMP_seqs %>% group_by(MAMP_Sequence) %>% summarise(n=n()))
csp22_occurance_by_genera <- as.data.frame(csp22_MAMP_seqs %>% group_by(MAMP_Sequence, Genera) %>% summarise(n=n()))

csp22_occurance_by_genera <- reshape(csp22_occurance_by_genera, idvar = "MAMP_Sequence", timevar = "Genera", direction = "wide")
csp22_occurance_by_genera[is.na(csp22_occurance_by_genera)] <- 0



ha <- rowAnnotation("MAMP \nOccurance" = anno_barplot(csp22_occurance$n, width = unit(2, "cm")),
                    "Streptomyces" = csp22_occurance_by_genera$n.Streptomyces,
                    "Rhodococcus" = csp22_occurance_by_genera$n.Rhodococcus,
                    "Clavibacter" = csp22_occurance_by_genera$n.Clavibacter,
                    "Leifsonia" = csp22_occurance_by_genera$n.Leifsonia,
                    "Curtobacterium" = csp22_occurance_by_genera$n.Curtobacterium,
                    "Rathayibacter" = csp22_occurance_by_genera$n.Rathayibacter,
                    
                    
                    "Xanthomonas" = csp22_occurance_by_genera$n.Xanthomonas,
                    "Agrobacterium" = csp22_occurance_by_genera$n.Agrobacterium,
                    "Ralstonia" = csp22_occurance_by_genera$n.Ralstonia,
                    "Pseudomonas" = csp22_occurance_by_genera$n.Pseudomonas
                   )

ComplexHeatmap::Heatmap(empty_csp22_matrix,
                        #row_km = 5,
                        #row_km_repeats = 1000,
                        #column_km = ,
  
                        #column_split = 4,
                        right_annotation = ha,

                        show_column_dend = F,
                        row_dend_width = unit(4.5,"cm"),                       
                        show_row_names = T,
                        show_column_names = F, 
                        row_names_gp = gpar(fontsize = 2.5),
                        
                        row_dend_reorder = T,
                        column_dend_reorder = T)



############################################################################ 
# repeat on subset of data -> MAMP has had to occur at teast 10 times
############################################################################ 


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
                                  csp22_occurance_by_genera %>% filter(n.Pseudomonas > 10)
                                  
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





###################################################################################################

#hold_cspD <- c("MQNGKVKWFNNEKGFGFIEVEGGDDVFVHFTAIEGDGYKSLEEGQEVSFEIVEGNRGPQASNVVKL")
#hold_cspB <- c("MLEGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTKEA")

test2 <- Biostrings::pairwiseAlignment(hold_cspD, All_target_by_annotation$seq[4], type = "global")
test2
Biostrings::pid(test2, type = "PID1")
