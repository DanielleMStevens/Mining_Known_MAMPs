######################################################################
# upload raw data 
######################################################################

# load all files
load_protein_whole_genome <- list.files(path = "./../../Whole_Genomes_v2/whole_genomes/", pattern = ".fna")


# will use datasettable to parse accession number for each genera to run ANI analysis on

for (i in 1:length(name_list)){
  accession_numbers_list <- list()
  subset_accesion_numbers <- subset(datasettable, datasettable$Genera == name_list[[1]])
  for (j in 1:nrow(subset_accesion_numbers)){
    which_genome <- grepl(subset_accesion_numbers$Assembly_Accession[1], load_protein_whole_genome)
    which_genome <- load_protein_whole_genome[which_genome]
    which_genome <- paste("./../../Whole_Genomes_v2/whole_genomes/", which_genome)
    accession_numbers_list[[j]] <- which_genome
  }
  
}