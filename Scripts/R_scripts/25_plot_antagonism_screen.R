

import_antagonist_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Figure6BD_Antagonism_Screen.xlsx", sheetName = "Summary")
colnames(import_antagonist_ros_data) <- c("Candidate Antagonist","Agonist","Water","Control","100 nM","200 nM","500 nM","1 uM")

import_antagonist_ros_data <- import_antagonist_ros_data[,c(1,2,4:8)]
import_antagonist_ros_data$Candidate_Antagonist_F = factor(import_antagonist_ros_data$`Candidate Antagonist`, levels=c("M-CP-43","M-CP-51","M-CP-59","M-CP-55","M-CP-53","Cm csp22-3","s-csp22-3"))
import_antagonist_ros_data <- reshape2::melt(import_antagonist_ros_data[,c(2:8)], id = c("Candidate_Antagonist_F","Agonist"))


(ggplot(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "con. csp22"), 
       aes(x = variable, y = `Candidate_Antagonist_F`)) +
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0, 100000), range = c(1,10), breaks = c(10000,20000,40000, 80000)) +
  xlab("") +
  ylab("Candidate Antagonist") +
  my_ggplot_theme +
    theme(axis.text.x = element_blank())) /
  
  ggplot(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "Cm csp22-1"), 
         aes(x = variable, y = `Candidate_Antagonist_F`)) +
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0, 100000), range = c(1,10), breaks = c(10000,20000,40000, 80000)) +
  xlab("\nAvg. Max RLUs") +
  ylab("Candidate Antagonist") +
  my_ggplot_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  patchwork::plot_layout(heights = c(3, 1))


test <- as.matrix(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "con. csp22")[c(3:8)])
rownames(test) <- unlist(as.list(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "con. csp22")[1]))

test2 <- as.matrix(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "Cm csp22-1")[c(3:8)])
rownames(test2) <- unlist(as.list(subset(import_antagonist_ros_data, import_antagonist_ros_data$Agonist == "Cm csp22-1")[1]))

(pheatmap(test,cluster_rows=FALSE, cluster_cols=FALSE, border_color = "white")) /
pheatmap(test2,cluster_rows=FALSE, cluster_cols=FALSE, border_color = "white")


# Export heatmap (4.5?)

Cm_csp22_3_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Antagonism_screen.xlsx", sheetName = "Cm csp22-3 | Cm csp22-1")
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
Cm_csp22_3_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Antagonism_screen.xlsx", sheetName = "con. elf18 | Cm csp22-1")
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
Cm_csp22_3_ros_data <- xlsx:::read.xlsx("./../../Analyses/MAMP_antagonism/Antagonism_screen.xlsx", sheetName = "s-csp22-3 | Cm csp22-1")
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
