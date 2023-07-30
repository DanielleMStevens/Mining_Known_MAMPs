#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


directory_path <- c('~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/')
type_fasta_files <- list.files(path = directory_path)


hold_type_csp_results <- data.frame("query_id" = character(0),
                                    "subject_id" = character(0),
                                    "perc_identity" = numeric(0),
                                    "num_ident_matches" = numeric(0),
                                    "alig_length"= numeric(0),
                                    "mismatches"= numeric(0),
                                    "gap_openings"= numeric(0),
                                    "n_gaps"= numeric(0),
                                    "pos_match"= numeric(0),
                                    "ppos"= numeric(0),
                                    "q_start"= numeric(0),
                                    "q_end"= numeric(0),
                                    "q_len"= numeric(0),
                                    "qcov"= numeric(0),
                                    "qcovhsp"= numeric(0),
                                    "s_start" = numeric(0),
                                    "s_end"= numeric(0),
                                    "s_len"= numeric(0),
                                    "evalue"= numeric(0),
                                    "bit_score" = numeric(0),
                                    "score_raw" = numeric(0),
                                    "Genus1" = character(0),
                                    "Type_v1" = character(0),
                                    "Genus2" = character(0),
                                    "Type_v2" = character(0))



for (i in 1:length(type_fasta_files)){
  for (j in 1:length(type_fasta_files)){
    genus <- str_split(type_fasta_files[[i]],"_")[[1]][1]
    type <- str_split(type_fasta_files[[i]],"_")[[1]][2]
    type <-str_split(type,"\\.")[[1]][1]
    
    genus2 <- str_split(type_fasta_files[[j]],"_")[[1]][1]
    type2 <- str_split(type_fasta_files[[j]],"_")[[1]][2]
    type2 <-str_split(type2,"\\.")[[1]][1]
    
    print(type_fasta_files[[i]])
    print(type_fasta_files[[j]])

    hold_df <- data.frame(orthologr::blast(query_file = paste0(directory_path, type_fasta_files[[i]], ""),
                                           subject_file = paste0(directory_path, type_fasta_files[[j]], ""),
                                           seq_type = 'protein',
                                           comp_cores = 8,
                                           eval = "1E0"))


    hold_df["Genus1"] <- rep(genus, nrow(hold_df))
    hold_df["Type_v1"] <- rep(type, nrow(hold_df))
    
    hold_df["Genus2"] <- rep(genus2, nrow(hold_df))
    hold_df["Type_v2"] <- rep(type2, nrow(hold_df))

    hold_type_csp_results <- rbind(hold_type_csp_results, hold_df)
    
  }
}

hold_type_csp_results$v1 <- paste(hold_type_csp_results$Genus1, hold_type_csp_results$Type_v1, sep = "_")
hold_type_csp_results$v2 <- paste(hold_type_csp_results$Genus2, hold_type_csp_results$Type_v2, sep = "_")


type_comp_sum <- data.frame(hold_type_csp_results %>% group_by(v1, v2) %>% summarise(m = mean(perc_identity)))


### reorganize table for plotting 
type_comp_Names <- sort(unique(as.character(unlist(type_comp_sum[2]))))
type_comp_matrix <- matrix(0, nrow = length(type_comp_Names), ncol = length(type_comp_Names),
                             dimnames = list(type_comp_Names, type_comp_Names))

type_comp_matrix[as.matrix(type_comp_sum[1:2])] <- type_comp_sum$m



Heatmap(type_comp_matrix,
        border = T,
        cluster_columns = T,
        cluster_rows = T,
        rect_gp = gpar(col = "white", lwd = 2))




#----------------------------------------------------------------------------------------------------------------------------------

# Clavibacter Typeing Comparison by Reciropcal Blast
Clavi_T1_T1 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Clavibacter_Type1.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Clavibacter_Type1.fasta', 
                                              seq_type = 'protein'))


Clavi_T2_T2 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Clavibacter_Type2.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Clavibacter_Type2.fasta', 
                                              seq_type = 'protein'))


Clavi_T1_T1 <- cbind(Clavi_T1_T1, "Comparison" = rep("Type 1 vs. Type 1", nrow(Clavi_T1_T1)))
Clavi_T2_T2 <- cbind(Clavi_T2_T2, "Comparison" = rep("Type 2 vs. Type 2", nrow(Clavi_T2_T2)))


Clavi_type_comp <- rbind(Clavi_T1_T1,Clavi_T2_T2)

ggplot(Clavi_type_comp, aes(x = Comparison, y = perc_identity)) +
  geom_quasirandom(method = "tukeyDense", size = 1, alpha = 0.5) +
  my_ggplot_theme +
  ylab("Percent Identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(20, 100)



# Leifsonia Typeing Comparison by Reciropcal Blast
Leif_T1_T1 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Leifsonia_Type1.fasta',
                                             subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Leifsonia_Type1.fasta', 
                                             seq_type = 'protein'))


Leif_T2_T2 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Leifsonia_Type2.fasta',
                                             subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Leifsonia_Type2.fasta', 
                                             seq_type = 'protein'))



Leif_T1_T1 <- cbind(Leif_T1_T1, "Comparison" = rep("Type 1 vs. Type 1", nrow(Leif_T1_T1)))
Leif_T2_T2 <- cbind(Leif_T2_T2, "Comparison" = rep("Type 2 vs. Type 2", nrow(Leif_T2_T2)))


Leif_type_comp <- rbind(Leif_T1_T1, Leif_T2_T2)


ggplot(Leif_type_comp, aes(x = Comparison, y = perc_identity)) +
  geom_quasirandom(method = "tukeyDense", size = 1, alpha = 0.5) +
  my_ggplot_theme +
  ylab("Percent Identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(20, 100)



# Rhodococcus Typeing Comparison by Reciropcal Blast
Rhodo_T1_T1 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type1.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type1.fasta', 
                                              seq_type = 'protein'))

Rhodo_T2_T2 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type2.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type2.fasta', 
                                              seq_type = 'protein'))

Rhodo_T3_T3 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type3.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type3.fasta', 
                                              seq_type = 'protein'))

Rhodo_T4_T4 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type4.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Rhodococcus_Type4.fasta', 
                                              seq_type = 'protein'))




Rhodo_T1_T1 <- cbind(Rhodo_T1_T1, "Comparison" = rep("Type 1 vs. Type 1", nrow(Rhodo_T1_T1)))
Rhodo_T2_T2 <- cbind(Rhodo_T2_T2, "Comparison" = rep("Type 2 vs. Type 2", nrow(Rhodo_T2_T2)))
Rhodo_T3_T3 <- cbind(Rhodo_T3_T3, "Comparison" = rep("Type 3 vs. Type 3", nrow(Rhodo_T3_T3)))
Rhodo_T4_T4 <- cbind(Rhodo_T4_T4, "Comparison" = rep("Type 4 vs. Type 4", nrow(Rhodo_T4_T4)))


Rhodo_type_comp <- rbind(Rhodo_T1_T1, Rhodo_T2_T2,
                         Rhodo_T3_T3, Rhodo_T4_T4)

ggplot(Rhodo_type_comp, aes(x = Comparison, y = perc_identity)) +
  geom_quasirandom(method = "tukeyDense", size = 1, alpha = 0.5) +
  my_ggplot_theme +
  ylab("Percent Identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(20, 100)




# Erwinina Typeing Comparison by Reciropcal Blast
Erwin_T1_T1 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Erwinia_Type1.fasta',
                                             subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Erwinia_Type1.fasta', 
                                             seq_type = 'protein'))


Erwin_T2_T2 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Erwinia_Type2.fasta',
                                             subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Erwinia_Type2.fasta', 
                                             seq_type = 'protein'))



Erwin_T1_T1 <- cbind(Erwin_T1_T1, "Comparison" = rep("Type 1 vs. Type 1", nrow(Erwin_T1_T1)))
Erwin_T2_T2 <- cbind(Erwin_T2_T2, "Comparison" = rep("Type 2 vs. Type 2", nrow(Erwin_T2_T2)))


Erwin_type_comp <- rbind(Erwin_T1_T1, Erwin_T2_T2)


ggplot(Erwin_type_comp, aes(x = Comparison, y = perc_identity)) +
  geom_quasirandom(method = "tukeyDense", size = 1, alpha = 0.5) +
  my_ggplot_theme +
  ylab("Percent Identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(20, 100)



# Xanth Typeing Comparison by Reciropcal Blast
Xanth_T1_T1 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type1.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type1.fasta', 
                                              seq_type = 'protein'))


Xanth_T2_T2 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type2.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type2.fasta', 
                                              seq_type = 'protein'))

Xanth_T4_T4 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type4.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type4.fasta', 
                                              seq_type = 'protein'))

Xanth_T5_T5 <- as.data.frame(orthologr::blast(query_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type5.fasta',
                                              subject_file = '~/Documents/Mining_MAMPs/Mining_Known_MAMPs/Analyses/Catagorizing_CSPs/CSP_types/Typing_fasta_files/Xanthomonas_Type5.fasta', 
                                              seq_type = 'protein'))


Xanth_T1_T1 <- cbind(Xanth_T1_T1, "Comparison" = rep("Type 1 vs. Type 1", nrow(Xanth_T1_T1)))
Xanth_T2_T2 <- cbind(Xanth_T2_T2, "Comparison" = rep("Type 2 vs. Type 2", nrow(Xanth_T2_T2)))
Xanth_T4_T4 <- cbind(Xanth_T4_T4, "Comparison" = rep("Type 4 vs. Type 4", nrow(Xanth_T4_T4)))
Xanth_T5_T5 <- cbind(Xanth_T5_T5, "Comparison" = rep("Type 5 vs. Type 5", nrow(Xanth_T5_T5)))


Xanth_type_comp <- rbind(Xanth_T1_T1, Xanth_T2_T2, Xanth_T4_T4, Xanth_T5_T5)


ggplot(Xanth_type_comp, aes(x = Comparison, y = perc_identity)) +
  geom_quasirandom(method = "tukeyDense", size = 1, alpha = 0.5) +
  my_ggplot_theme +
  ylab("Percent Identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(20, 100)




















