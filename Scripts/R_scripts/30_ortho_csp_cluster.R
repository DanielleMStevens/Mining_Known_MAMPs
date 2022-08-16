# filtering and grouping csp cluster from mmseq2


mmseq_clusterin_linclust <- read.table(file = file.choose(), header= F, stringsAsFactors = F)
mmseq_clusterin_cluster <- read.table(file = file.choose(), header= F, stringsAsFactors = F)



colnames(mmseq_clusterin_cluster) <- c("group_reps","group_members")
mmseq_clusterin_cluster[,3:8] <- stringr::str_split_fixed(mmseq_clusterin_cluster$group_members, "\\|",6)
colnames(mmseq_clusterin_cluster) <- c("group_reps","group_members", "WP_tag","MAMP","ignore","Genus","FileName","ignore2")


for (i in 1:length(unique(mmseq_clusterin_cluster$group_reps))){
  test_rep_group <- subset(mmseq_clusterin_cluster, mmseq_clusterin_cluster$group_reps == unique(mmseq_clusterin_cluster$group_reps)[i])
  All_target_by_annotation_test <- All_target_by_annotation[All_target_by_annotation$Filename %in% test_rep_group$FileName,]
  All_target_by_annotation_test <- All_target_by_annotation_test[All_target_by_annotation_test$Protein_Name %in% test_rep_group$WP_tag,]                                              
  
  
  All_target_by_annotation_test <- formate2fasta(All_target_by_annotation_test$Protein_Name, paste0("csp22_type", i), All_target_by_annotation_test$Genera, 
                                                 All_target_by_annotation_test$Filename, All_target_by_annotation_test$seq)
  writeFasta(All_target_by_annotation_test, paste0("./../Protein_alignments_and_trees/csp_type_trees/","csp_type", i ,".fasta"))
}
