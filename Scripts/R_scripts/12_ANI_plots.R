#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 12/20/2020
# Script Purpose: Processing and Plotting ANI values
# Inputs Necessary: ANI_comp_name_sheet.txt, output from fast ANI analysis
# Outputs: plots all by all comparison of genomes, organizes by clustering of ANI values, and labels based on species 
#-----------------------------------------------------------------------------------------------

# moved values to new variable such that orginal dataframe is not altered just in case

hold_ANI_values <- as.data.frame(ANI_values_df)
hold_ANI_values <- data.frame(lapply(hold_ANI_values, unlist))
rownames(hold_ANI_values) <- NULL


#########################################################################
# remove genomes from genomes_to_check that are clonal
#########################################################################

hold_ANI_values <- hold_ANI_values[!hold_ANI_values$Genome1 %in% Genomes_to_check$Genome1,]
hold_ANI_values <- hold_ANI_values[!hold_ANI_values$Genome2 %in% Genomes_to_check$Genome1,]


#########################################################################
# reorganize table for plotting
#########################################################################

myNames <- sort(unique(as.character(unlist(hold_ANI_values[1:2]))))

empty_ANI_matrix <- matrix(0, 
                           nrow = length(unique(hold_ANI_values$Genome1)), 
                           ncol = length(unique(hold_ANI_values$Genome2)),
                           dimnames = list(myNames, myNames))

# fill in upper triangle
empty_ANI_matrix[as.matrix(hold_ANI_values[c(1,2)])] <- hold_ANI_values$ANI_val


#########################################################################
# define colors + plot heatmap of ANI values using complexHeatmap package
#########################################################################


set_colors_df <- datasettable[datasettable$Filename %in% unique(hold_ANI_values$Genome1),c(6,1)]
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
                                       
                                       use_raster = TRUE, raster_quality = 5,
                                       
                                       heatmap_legend_param = list(
                                         title = "ANI Value",
                                         legend_height = unit(3, "cm"),
                                         legend_width = unit(2.5, "cm"),
                                         border = "black",
                                         direction = "vertical"
                                       )
)


pdf("./../Figures/Supplemental_Figure_3/ANI_heatmap.pdf", width = 12, height = 12)
ANI_heatmap
dev.off()



#########################################################################
# subset each group by Genera and plot values
#########################################################################


copy_ANI <- hold_ANI_values
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

# number of strains left in study
n_number <- as.data.frame(set_colors_df %>% group_by(Genera) %>% summarise(n=n()))

#order of ggplot to match ANI heatmap
level_order <- c('Xanthomonas', 'Pseudomonas', 'Agrobacterium', 'Ralstonia', 'Pectobacterium', 'Erwinia', 
                 'Curtobacterium', 'Dickeya', 'Streptomyces', 'Rhodococcus', 'Clavibacter', 'Leifsonia', 'Rathayibacter')


ANI_plot_ggplot <- ggplot(copy_ANI, 
                          aes(x = factor(Genera, level = level_order), 
                              y = ANI_val,fill = Genera, color = Genera)) +
  geom_violin(alpha = 0.4, scale = "width", trim = T) +
  geom_boxplot(color = "black", fill = "white", outlier.alpha = 0, width = 0.05) +
  scale_y_continuous(breaks = c(75,80,85,90,95,100), limits = c(75,110)) +
  xlab("Genera\n") +
  ylab("ANI Value") +
  scale_color_manual("Genera", values = Genera_colors) +
  scale_fill_manual("Genera", values = Genera_colors) +
  my_ggplot_theme +
  theme(axis.text.x = element_text(color ="black", angle = 45, hjust = 1),
        legend.position = "none") +
  geom_text(data = n_number, 
            aes(x = Genera, y = 105, label = n), color = "black", size = 4) 



ggsave(ANI_plot_ggplot, filename = "./../Figures/Supplemental_Figure_3/ANI_ggplot.pdf", device = cairo_pdf, width = 7.5, height = 3, units = "in")


