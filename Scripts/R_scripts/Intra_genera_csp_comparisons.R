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



