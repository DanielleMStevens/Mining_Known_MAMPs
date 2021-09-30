######################################################################
# upload raw data 
######################################################################

# load all files
load_protein_whole_genome <- list.files(path = "./../../Whole_Genomes_v2/whole_genomes/", pattern = ".fna")


# will use datasettable to parse accession number for each genera to run ANI analysis on

for (i in 1:length(name_list)){
  accession_numbers_list <- list()
  subset_accesion_numbers <- subset(datasettable, datasettable$Genera == name_list[[i]])
  for (j in 1:nrow(subset_accesion_numbers)){
    which_genome <- grepl(subset_accesion_numbers$Assembly_Accession[j], load_protein_whole_genome)
    which_genome <- load_protein_whole_genome[which_genome]
    which_genome <- paste("/home/danimstevens/Documents/Mining_MAMPs/Whole_Genomes_v2/whole_genomes/", which_genome, sep = "")
    accession_numbers_list[[j]] <- which_genome
  }
  writeLines(unlist(accession_numbers_list), paste("/home/danimstevens/Documents/Mining_MAMPs/Mining_Known_MAMPs/ANI_analysis/ANI_Analysis_", name_list[i], ".txt", sep = ""))
}
