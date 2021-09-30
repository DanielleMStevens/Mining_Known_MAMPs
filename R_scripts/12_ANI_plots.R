
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
                                       show_row_names = T,
                                       row_names_gp = gpar(fontsize = 2.5),
                                       
                                       width = unit(20, "cm"),
                                       height = unit(20, "cm"),
                                       
                                       # left_annotation = ha,
                                       #  top_annotation = column_ha,
                                       
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
        legend.position = "none") 

coord_flip()


ANI_plot_ggplot

#########################################################################
# subset each group by Genera and plot values
#########################################################################

library(cowplot)
cowplot::plot_grid(ANI_heatmap, ANI_plot_ggplot, ncol = 1, align = "v")


cowplot::plot_grid()