#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


##############################################
# importing MAMP Screen data
##############################################

MAMP_to_test_elf18 <- as.data.frame(read_xlsx("./../ROS_Screen/MAMPs_Tested.xlsx", sheet = 1))
#MAMP_to_test_elf18 <- MAMP_to_test[73:97,]

elf18_sim_score_to_test <- data.frame("MAMP_v1" = character(0), "MAMP_v2" = character(0), "Percent_sim" = numeric(0))
for (i in 1:nrow(MAMP_to_test_elf18)){
  for (j in 2:nrow(MAMP_to_test_elf18)){
    elf18_sim_score_to_test <- rbind(elf18_sim_score_to_test, data.frame("MAMP_v1" = MAMP_to_test_elf18$Sequence[[i]],
                                                         "MAMP_v2" = MAMP_to_test_elf18$Sequence[[j]],
                                                         "Percent_sim" = Biostrings::pid(Biostrings::pairwiseAlignment(MAMP_to_test_elf18$Sequence[[i]], 
                                                                                                                       MAMP_to_test_elf18$Sequence[[j]], 
                                                                                                                       type = "global", 
                                                                                                                       substitutionMatrix = BLOSUM62), 
                                                                                         type = "PID1")))
  }
} 

# fill in upper triangle
elf18_Names <- sort(unique(as.character(unlist(elf18_sim_score_to_test[1:2]))))
empty_elf18_matrix <- matrix(0, nrow = length(elf18_Names), ncol = length(elf18_Names),dimnames = list(elf18_Names, elf18_Names))
empty_elf18_matrix[as.matrix(elf18_sim_score_to_test[c(1,2)])] <- elf18_sim_score_to_test$Percent_sim


row_dend = as.dendrogram(hclust(dist(empty_elf18_matrix)))
# Rectangular lines
#plot(as.phylo(row_dend), type = "fan", cex = 0.6, label.offset = 0.5)


hc <- hclust(dist(empty_elf18_matrix), "ave")
hcdata <- dendro_data(hc, type = "rectangle")
ggplot() +
  geom_segment(data = segment(hcdata), 
               aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  geom_text(data = label(hcdata), 
            aes(x = x, y = y, label = label, hjust = 0), 
            size = 3
  ) +
  coord_flip() +
  scale_y_reverse(expand = c(0.8, 0)) + 
  my_ggplot_theme +
  theme(rect = element_blank(),
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.title.x = element_blank(),
       axis.title.y = element_blank())


