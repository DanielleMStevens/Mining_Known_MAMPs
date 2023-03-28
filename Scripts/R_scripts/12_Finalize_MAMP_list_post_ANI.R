#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


######################################################################
# fix the double domain hits in CSPs
######################################################################

#essentially (and became apparent in Agro vitus S4), some CSPS have two CSP domains and rather than have each listed, it
#is currently the same on twice. The script below should fix that

which_CSPS_are_double_domains <- subset(filtered_hold_MAMP_seqs, filtered_hold_MAMP_seqs$MAMP_Hit == "csp22_consensus")
which_CSPS_are_double_domains <- as.data.table(which_CSPS_are_double_domains %>% group_by(File_Name, Protein_Name) %>% summarise(n=n()))
which_CSPS_are_double_domains <- subset(which_CSPS_are_double_domains, which_CSPS_are_double_domains$n > 1)


hold_corrected_CSP_doublets <- data.frame("File_Name" = character(0), "Protein_Name" = character(0), "First_version"= character(0), "Second_version" = character(0))
for (i in 1:nrow(which_CSPS_are_double_domains)) {
  subset_the_genome <- subset(All_target_by_annotation, All_target_by_annotation$Filename == which_CSPS_are_double_domains$File_Name[i])
  subset_the_gene <- subset(subset_the_genome, subset_the_genome$Protein_Name == which_CSPS_are_double_domains$Protein_Name[i])
  second_hit <- Biostrings::matchPattern(load_reference_MAMPs_fasta$seq[1], subset_the_gene$seq, max.mismatch = 15)
  hold_corrected_CSP_doublets <- rbind(hold_corrected_CSP_doublets,
                                       data.frame("File_Name" = subset_the_gene$Filename,
                                                  "Protein_Name" = subset_the_gene$Protein_Name,
                                                  "First_version" = substr(subset_the_gene$seq, second_hit@ranges@start[1], (second_hit@ranges@start[1]+ second_hit@ranges@width) - 1),
                                                  "Second_version" = substr(subset_the_gene$seq, second_hit@ranges@start[2], (second_hit@ranges@start[2]+second_hit@ranges@width) - 1)
                                       ))
}
rm(subset_the_gene)
rm(subset_the_genome)
rm(second_hit)

for (i in 1:nrow(hold_corrected_CSP_doublets)){
  for (j in 1:nrow(filtered_hold_MAMP_seqs)){
    if (hold_corrected_CSP_doublets$File_Name[i] == filtered_hold_MAMP_seqs$File_Name[j]){
      if (hold_corrected_CSP_doublets$Protein_Name[i] == filtered_hold_MAMP_seqs$Protein_Name[j]){
        if(filtered_hold_MAMP_seqs[j,4] == hold_corrected_CSP_doublets$First_version[i]){
          filtered_hold_MAMP_seqs[j+1,4] <- hold_corrected_CSP_doublets$Second_version[i]
          filtered_hold_MAMP_seqs[j+1,3] <- pid(pairwiseAlignment(load_reference_MAMPs_fasta$seq[1], hold_corrected_CSP_doublets$Second_version[i]), type = "PID1")
        }
        #print(filtered_hold_MAMP_seqs[j,])
      }
    }
  }  
}

rm(hold_corrected_CSP_doublets)
rm(which_CSPS_are_double_domains)


# fix naming issues to be congruent for tree building (run 6 times for remove all :,;# symbols)
filtered_hold_MAMP_seqs$File_Name <- stringr::str_replace(filtered_hold_MAMP_seqs$File_Name, '\\:', '\\_')
filtered_hold_MAMP_seqs$File_Name <- stringr::str_replace(filtered_hold_MAMP_seqs$File_Name, '\\:', '\\_')
filtered_hold_MAMP_seqs$File_Name <- stringr::str_replace(filtered_hold_MAMP_seqs$File_Name, '\\:', '\\_')
filtered_hold_MAMP_seqs$File_Name <- stringr::str_replace(filtered_hold_MAMP_seqs$File_Name, '\\,', '\\_')
filtered_hold_MAMP_seqs$File_Name <- stringr::str_replace(filtered_hold_MAMP_seqs$File_Name, '\\;', '\\_')
filtered_hold_MAMP_seqs$File_Name <- stringr::str_replace(filtered_hold_MAMP_seqs$File_Name, '\\#', '\\_')


######################################################################
# calculate copy number
######################################################################


hold_copy_number <- as.data.frame(filtered_hold_MAMP_seqs %>% group_by(MAMP_Hit, File_Name, Gram, Genera) %>% summarise(n=n()))


######################################################################
# create the final dataframe into an html file so others can 
######################################################################

DT::datatable(filtered_hold_MAMP_seqs, 
              colnames = c('WP Tag', 'MAMP Hit', 'Percent Identity', 
                           'MAMP Sequence', 'Genera', 'Strain', 
                           'File Name', 'Gram Type'),
              rownames = FALSE,
              filter = 'top',
              options = list(
                autoWidth = TRUE,
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
      "}"))
)



