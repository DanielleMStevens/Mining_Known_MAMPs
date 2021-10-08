#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 12/20/2020
# Script Purpose: Processing and Plotting ANI values
# Inputs Necessary: ANI_comp_name_sheet.txt, output from fast ANI analysis
# Outputs: plots all by all comparison of genomes, organizes by clustering of ANI values, and labels based on species 
#-----------------------------------------------------------------------------------------------

# moved values to new variable such that orginal dataframe is not altered just in case

hold_ANI_values <- ANI_values_df


#########################################################################
# remove genomes from genomes_to_check that are clonal
#########################################################################

test <- hold_ANI_values[!hold_ANI_values$Genome1 %in% Genomes_to_check$Genome1,]
test <- test[!test$Genome2 %in% Genomes_to_check$Genome1,]



#########################################################################
# reorganize table for plotting
#########################################################################

myNames <- sort(unique(as.character(unlist(test[1:2]))))

empty_ANI_matrix <- matrix(0, nrow = length(unique(test$Genome1)), 
                           ncol = length(unique(test$Genome2)),
                           dimnames = list(myNames, myNames))

# fill in upper triangle
empty_ANI_matrix[as.matrix(test[c(1,2)])] <- test$ANI_val
# fill in the lower triangle
empty_ANI_matrix[as.matrix(test[c(2,1)])] <- test$ANI_val

#

#########################################################################
# define colors + plot heatmap of ANI values using complexHeatmap package
#########################################################################


set_colors_df <- datasettable[datasettable$Filename %in% unique(test$Genome1),c(6,1)]
set_colors_df <- set_colors_df[!duplicated(set_colors_df$Filename),]
set_colors_df <- set_colors_df[match(rownames(empty_ANI_matrix), set_colors_df$Filename),]

########################################################################
# reorganize table for plotting
#########################################################################


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
                                       #row_names_gp = gpar(fontsize = 2.5),
                                       
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


#for (i in nrow(copy_ANI):1){
#  if (strsplit(copy_ANI$Genome1[i],"_")[[1]][1] == strsplit(copy_ANI$Genome2[i],"_")[[1]][1]){
#    next
#  }
# if (strsplit(copy_ANI$Genome1[i],"_")[[1]][1] != strsplit(copy_ANI$Genome2[i],"_")[[1]][1]){
#    copy_ANI[i,] <- NA
#  }
#}


#copy_ANI <- na.omit(copy_ANI)

copy_ANI <- test
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

n_number <- as.data.frame(copy_ANI %>% group_by(Genera) %>% summarise(n=n()))

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
        legend.position = "none") 

coord_flip()


ANI_plot_ggplot

