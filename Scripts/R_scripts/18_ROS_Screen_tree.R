#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: After screeening all the elf18 and csp22 peptides, I wanted to guide each ROS screen based onthe dervied peptide data itself using basic phylogeic trees 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


##############################################
# importing phylogenetic tree - elf18 ROS variants
##############################################

elf18_ROS_tree <- read.tree("./../../Analyses/ROS_Screen/eptiope_tree_for_ROS/elf18_variants_tree.tre")
elf18_ROS_tree <- phangorn::midpoint(elf18_ROS_tree, node.labels='label')

ggtree::ggtree(elf18_ROS_tree,  ladderize = T, size = 0.2, linetype = 1, branch.length = "none") + 
  geom_tiplab(aes(label = label), offset = 7, hjust = 1, size = 2.5)


# image was exported as pdf at 1.5 x 3.3 inch size.

csp22_ROS_tree <- read.tree("./../../Analyses/ROS_Screen/eptiope_tree_for_ROS/csp22_variants_tree.tre")
csp22_ROS_tree <- phangorn::midpoint(csp22_ROS_tree, node.labels='label')

ggtree::ggtree(csp22_ROS_tree,  ladderize = T, size = 0.5, linetype = 1,  branch.length = "none") + 
  geom_tiplab(aes(label = label), offset = 14, hjust = 1, size = 2.5)


# image was exported as pdf at 2 x 7.5 inch size.


##############################################
# importing ROS Screen outputs as data tables to import onto phylogenetic tree
##############################################


import_ros_screen_outcome <- xlsx:::read.xlsx("./../../Analyses/ROS_Screen/ROS_Screen_data.xlsx", sheetName = "Epitope_Immunogenecity_Outcomes")
elf18_ros_screen_outcome <- import_ros_screen_outcome[69:94,]
csp22_ros_screen_outcome <- rbind(import_ros_screen_outcome[1:68,], import_ros_screen_outcome[95:97,])



##############################################
# Assessing how similarity equates to immunogenicity
##############################################

#elf18_ros_screen_outcome_ggplot <- elf18_ros_screen_outcome
#elf18_ros_screen_outcome_ggplot["Percent_Identity"] <- numeric(0)
#for (i in 1:nrow(elf18_ros_screen_outcome_ggplot)){
#  elf18_ros_screen_outcome_ggplot$Percent_Identity[i] <- unique(subset(filtered_hold_MAMP_seqs, 
#                                                                       filtered_hold_MAMP_seqs$MAMP_Sequence == elf18_ros_screen_outcome_ggplot$Sequence[i])[3])[[1]]
#}

#ggplot(elf18_ros_screen_outcome_ggplot, aes(x = 1, y = Percent_Identity, color = Outcome)) +
#  geom_jitter(width = 0.2) +
#  my_ggplot_theme


#csp22_ros_screen_outcome_ggplot <- csp22_ros_screen_outcome
#csp22_ros_screen_outcome_ggplot["Percent_Identity"] <- numeric(0)
#for (i in 1:nrow(csp22_ros_screen_outcome_ggplot)){
#  csp22_ros_screen_outcome_ggplot$Percent_Identity[i] <- unique(subset(filtered_hold_MAMP_seqs, 
#                                                                       filtered_hold_MAMP_seqs$MAMP_Sequence == csp22_ros_screen_outcome_ggplot$Sequence[i])[3])[[1]]
#}
#csp22_ros_screen_outcome_ggplot["MAMP_type"] <- c("csp22_variant")

#ggplot(csp22_ros_screen_outcome_ggplot, aes(x = factor(Outcome, levels = c("Immunogenic","Slightly Immunogenic","Non-Immunogenic","Not Tested")), 
#                                            y = Percent_Identity, color = Outcome)) +
#  geom_quasirandom(method = "tukeyDense", size = 0.5) +
#  geom_boxplot(alpha = 0.7, outlier.fill = NA, width = 0.2) +
#  scale_color_manual("", values = Immunogenicity_colors) +
#  my_ggplot_theme +
#  xlab("") +
#  ylab("Percent AA Similarity") +
#  ylim(0,100) +
#  coord_flip() +
#  theme(legend.position = "none", axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7),
#        axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))




################ old code. saving for future reference ########################

#hc <- hclust(dist(empty_elf18_matrix), "ave")
#hcdata <- dendro_data(hc, type = "rectangle")
#ggplot() +
#  geom_segment(data = segment(hcdata), 
#               aes(x = x, y = y, xend = xend, yend = yend)
#  ) +
#  geom_text(data = label(hcdata), 
#            aes(x = x, y = y, label = label, hjust = 0), 
#            size = 3
#  ) +
#  coord_flip() +
#  scale_y_reverse(expand = c(0.8, 0)) + 
#  my_ggplot_theme +
#  theme(rect = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank())


