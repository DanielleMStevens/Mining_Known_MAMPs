#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 12/20/2020
# Script Purpose: Processing and Plotting ANI values
# Inputs Necessary: ANI_comp_name_sheet.txt, output from fast ANI analysis
# Outputs: plots all by all comparison of genomes, organizes by clustering of ANI values, and labels based on species 
#-----------------------------------------------------------------------------------------------

#########################################################################
# load libraries
#########################################################################

library(ComplexHeatmap)
library(reshape2)
library(stringr)
library(cluster)

#########################################################################
# point to folder with DNA contig files, need to all be formated by Species_Strain.fasta (ex. CC_PF008.fasta)
#########################################################################

# defines a name list to consistently order genera
name_list <- c('Clavibacter','Leifsonia','Rathayibacter','Curtobacterium','Rhodococcus','Streptomyces',
               'Agrobacterium','Ralstonia','Xanthomonas','Pseudomonas')




if(exists("datasettable") == FALSE){
  source("./tree_tip_data.R")
}

########################################################################
# input files 
#########################################################################

matrix_table <- read.table(file = file.choose()) # select the output of the ANI analysis (file named ANI_analysis)


#########################################################################
# update table variables
#########################################################################


ANI_output <- matrix_table
colnames(ANI_output) <- c("Genome1","Genome2","ANI_val","val2","val3")
ANI_output$Genome1 <- as.character(ANI_output$Genome1)
ANI_output$Genome2 <- as.character(ANI_output$Genome2)
ANI_output$val2 <- NULL
ANI_output$val3 <- NULL

#########################################################################
# create a filter function, which will help remove information (aka the path)
# from the name of each strain such that it can be plotted with just its strain name
#########################################################################


filter_title <- function(which_table, which_phrase){
  for (i in 1:nrow(which_table)){
    which_table[i,1] <- stringr::str_replace(which_table[i,1], which_phrase, "")
    which_table[i,2] <- stringr::str_replace(which_table[i,2], which_phrase, "")
  }
  return(which_table) 
}


ANI_output <- filter_title(ANI_output,"/home/danimstevens/Documents/Mining_MAMPs/Whole_Genomes/whole_genome_fasta/")
ANI_output <- filter_title(ANI_output,".fna")


#########################################################################
# swap out accession number for the file name so we can map on genera metasdata on it
#########################################################################


 for (i in 1:nrow(ANI_output)){
   hold_genome1 <- paste(strsplit(ANI_output$Genome1[i], "_")[[1]][1],
                         strsplit(ANI_output$Genome1[i], "_")[[1]][2],
                         sep = "_")
   hold_genome2 <- paste(strsplit(ANI_output$Genome2[i], "_")[[1]][1],
                         strsplit(ANI_output$Genome2[i], "_")[[1]][2],
                         sep = "_")
   ANI_output$Genome1[i] <- datasettable[datasettable$Assembly_Accession %in% hold_genome1,7]
   ANI_output$Genome2[i] <- datasettable[datasettable$Assembly_Accession %in% hold_genome2,7]
   
 }


#########################################################################
# anything below 80 is equal to 0
#########################################################################


for (i in 1:nrow(ANI_output)){
  if (strsplit(ANI_output$Genome1[i], "_")[[1]][1] != strsplit(ANI_output$Genome2[i], "_")[[1]][1]){
    ANI_output$ANI_val[i] <- 0
  }
}


#########################################################################
# reorganize table for plotting
#########################################################################

myNames <- sort(unique(as.character(unlist(ANI_output[1:2]))))

empty_ANI_matrix <- matrix(0, nrow = length(unique(ANI_output$Genome1)), 
                           ncol = length(unique(ANI_output$Genome2)),
                           dimnames = list(myNames, myNames))

# fill in upper triangle
empty_ANI_matrix[as.matrix(ANI_output[c(1,2)])] <- ANI_output$ANI_val
# fill in the lower triangle
empty_ANI_matrix[as.matrix(ANI_output[c(2,1)])] <- ANI_output$ANI_val

#########################################################################
# reorganize table for plotting
#########################################################################



#########################################################################
# define colors + plot heatmap of ANI values using complexHeatmap package
#########################################################################


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

set_colors_df <- datasettable[datasettable$Filename %in% unique(ANI_output$Genome1),c(7,1)]
set_colors_df <- set_colors_df[!duplicated(set_colors_df$Filename),]
set_colors_df <- set_colors_df[match(rownames(empty_ANI_matrix), set_colors_df$Filename),]


ha <- rowAnnotation(Genera = set_colors_df$Genera,
                    border = T,
                    col = list(Genera = color_species))

column_ha <- HeatmapAnnotation(Genera2 = set_colors_df$Genera,
                               border = T,
                               show_legend = c("Genera2" = FALSE),
                               show_annotation_name = c("Genera2" = FALSE),
                               col = list(Genera2 = color_species))


row_dend = as.dendrogram(hclust(dist(empty_ANI_matrix)))


ANI_heatmap <- ComplexHeatmap::Heatmap(empty_ANI_matrix,
                                 #      cluster_rows = diana,
                              #cluster_rows = agnes,
                              #cluster_columns = agnes,
                              cluster_rows = row_dend,
                              cluster_columns = row_dend,
                              row_dend_reorder = T,
                              column_dend_reorder = T,
                              border = T, 
                              col = viridis::magma(nrow(empty_ANI_matrix)),
                              #row_km = 10,
                        
                              row_dend_width = unit(4.5,"cm"),
                              show_column_dend = F,
                              show_column_names = F, 
                              show_row_names = F,
                              width = unit(20, "cm"),
                              height = unit(20, "cm"),
                              
                              left_annotation = ha,
                              top_annotation = column_ha,
                              
                              heatmap_legend_param = list(
                                title = "ANI Value",
                                legend_height = unit(3, "cm"),
                                legend_width = unit(2.5, "cm"),
                                border = "black",
                                direction = "vertical"
                              )
)

ANI_heatmap

#########################################################################
# subset each group by Genera and plot values
#########################################################################


copy_ANI <- ANI_output
for (i in nrow(copy_ANI):1){
  if (strsplit(copy_ANI$Genome1[i],"_")[[1]][1] == strsplit(copy_ANI$Genome2[i],"_")[[1]][1]){
    next
  }
  if (strsplit(copy_ANI$Genome1[i],"_")[[1]][1] != strsplit(copy_ANI$Genome2[i],"_")[[1]][1]){
    copy_ANI[i,] <- NA
  }
}


copy_ANI <- na.omit(copy_ANI)


hold_genera_name <- list()
for (i in 1:nrow(copy_ANI)){
  hold_genera_name[[i]] <- strsplit(copy_ANI$Genome1[i],"_")[[1]][1]
}

copy_ANI <- cbind(copy_ANI, unlist(hold_genera_name))
colnames(copy_ANI) <- c("Genome1","Genome2","ANI_val","Genera")


# remove comparisons against itself (is. artifical 100%'s)
for (i in 1:nrow(copy_ANI)){
  if (copy_ANI$Genome1[i] == copy_ANI$Genome2[i]){
    print(i)
    copy_ANI[i,] <- NA
  }
}


copy_ANI <- na.omit(copy_ANI)


ANI_plot_ggplot <- ggplot(copy_ANI, aes(x = factor(Genera, level = name_list), y = ANI_val,
                    fill = Genera, color = Genera)) +
  geom_violin(alpha = 0.4, scale = "width", trim = F) +
  #geom_jitter() +
  ylim(75, 100) +
  xlab("Genera\n") +
  ylab("ANI Value\n") +
  scale_color_manual("Genera", values = Genera_colors) +
  scale_fill_manual("Genera", values = Genera_colors) +
  my_ggplot_theme +
  theme(axis.text.x = element_text(color ="black"),
        legend.position = "none") +
  coord_flip()


ANI_plot_ggplot

xs#########################################################################
# subset each group by Genera and plot values
#########################################################################

library(cowplot)
cowplot::plot_grid(ANI_heatmap, ANI_plot_ggplot, ncol = 1, align = "v")


cowplot::plot_grid()