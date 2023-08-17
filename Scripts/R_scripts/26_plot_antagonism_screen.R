#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 8/08/2023
# Script Purpose: Plot Antagonism Screen as a Heatmap
# Inputs: Figure6BD_Antagonism_Screen.xlsx 
# Outputs: Heatmaps for Antagonism Screen
#-----------------------------------------------------------------------------------------------

##############################################
# Import ROS data for antagonism screen
##############################################

import_antagonist_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Figure6BD_Antagonism_Screen.xlsx", sheetName = "Summary")
colnames(import_antagonist_ros_data) <- c("Candidate Antagonist","Agonist","Water","Control","100 nM","200 nM","500 nM","1 uM")


# plot those screened via consensus csp22 -------------- Figure6B ------------------
test <- as.matrix(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "con. csp22")[c(3:8)])
rownames(test) <- unlist(as.list(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "con. csp22")[1]))

pheatmap(test,cluster_rows=FALSE, cluster_cols=FALSE, border_color = "white")

# plot those screened via Cm csp22-1 -------------- Figure6D ------------------
test2 <- as.matrix(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "Cm csp22-1")[c(3:8)])
rownames(test2) <- unlist(as.list(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "Cm csp22-1")[1]))


pheatmap(test2,cluster_rows=FALSE, cluster_cols=FALSE, border_color = "white")


# Export Individual Replicates -------------- Figure6E ------------------

Cm_csp22_3_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Figure6BD_Antagonism_Screen.xlsx", sheetName = "Cm csp22-3 | Cm csp22-1")
colnames(Cm_csp22_3_ros_data) <- c("Replicate","Water","Control","100 nM","200 nM","500 nM","1 uM")

Cm_csp22_3_ros_data <- reshape2::melt(Cm_csp22_3_ros_data[1:16,])

(ggplot(subset(Cm_csp22_3_ros_data,Replicate == "04_07_2022_ROS_M06_all_conc_vs_Cm1_Max") , aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_point()) +

(ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_03_2022_Cm3_anti_Cm1_rep1") , aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_point())

ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_03_2022_Cm3_anti_Cm1_rep2") , aes(x = variable, y = value)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_point() +
  my_ggplot_theme +
  ylab("Max RLUs") +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# look at con. elf18
Cm_csp22_3_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Figure6BD_Antagonism_Screen.xlsx", sheetName = "con. elf18 | Cm csp22-1")
colnames(Cm_csp22_3_ros_data) <- c("Replicate","Water","Control","100 nM","200 nM","500 nM","1 uM")

Cm_csp22_3_ros_data <- reshape2::melt(Cm_csp22_3_ros_data[1:16,])

(ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_04_2022_elf18_anti_Cm1_rep1") , aes(x = variable, y = value)) +
    geom_boxplot() +
    geom_point()) 
  
  ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_15_2022_elf18_anti_Cm1") , aes(x = variable, y = value)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_point() +
  my_ggplot_theme +
  ylab("Max RLUs") +
  xlab("") + 
  ylim(0,120000) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


  
  ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_09_2022_elf18_anti_Cm1") , aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_point()



# look at s-csp22-3
Cm_csp22_3_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Figure6BD_Antagonism_Screen.xlsx", sheetName = "s-csp22-3 | Cm csp22-1")
colnames(Cm_csp22_3_ros_data) <- c("Replicate","Water","Control","100 nM","200 nM","500 nM","1 uM")

Cm_csp22_3_ros_data <- reshape2::melt(Cm_csp22_3_ros_data[1:12,])

ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_03_2022_scram_anti_Cm1_rep1") , aes(x = variable, y = value)) +
    geom_boxplot(fill = "grey", color = "black") +
    geom_point() +
    my_ggplot_theme +
    ylab("Max RLUs") +
    xlab("") + 
    ylim(0,125000) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  (ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_03_2022_scram_anti_Cm1_rep2") , aes(x = variable, y = value)) +
     geom_boxplot() +
     geom_point()) +
  
  ggplot(subset(Cm_csp22_3_ros_data,Replicate == "11_09_2022_scram_anti_Cm1") , aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_point()




