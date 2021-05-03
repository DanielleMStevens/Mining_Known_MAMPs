#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------



######################################################################
# set path to data
######################################################################

#setwd to where repo was cloned and maintained
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
#getSrcDirectory(function(x) {x})



##############################################
# Load colors - Set tip labels to match other figures
##############################################

#make sure to set path to the same place where the figure 
source("./package_dependencies.R")

source("./figure_colors.R")

source("./Theme_ggplot.R")

source("./CommonFunctions.R")

######################################################################
# upload raw data 
######################################################################

# load all files
load_protein_fasta_files <- list.files(path = "./../../Whole_Genomes/protein_fasta/", pattern = ".faa")
load_MAMP_blast_result_files <- list.files(path = './../MAMP_blast_outputs/', pattern = ".faa.txt")

# filter
load_protein_fasta_files <- load_protein_fasta_files[!load_protein_fasta_files %in% load_MAMP_blast_result_files]



#fix path so can load into process script
for (i in 1:length(load_protein_fasta_files)){
  load_protein_fasta_files[i] <- paste("./../../Whole_Genomes/protein_fasta/", load_protein_fasta_files[i], sep="")
  load_MAMP_blast_result_files[i] <- paste("./../MAMP_blast_outputs/", load_MAMP_blast_result_files[i], sep="")
}


# load reference file to determine percent identity of mined MAMPs
load_reference_MAMPs_fasta <- Biostrings::readAAStringSet(filepath = "./../MAMP_database/MAMP_elicitor_list.fasta") #may come with warning messgae..ignore for now
load_reference_MAMPs_fasta <- dss2df(load_reference_MAMPs_fasta)


##############################################
# load data which allows us to change accession name to strain name and add in genus info
##############################################


tip_data <- file.choose() #choose Mining_for_known_MAMPs_genome_accession_info.xlsx
datasettable <- readxl::read_xlsx(tip_data, col_names = T)
datasettable <- as.data.frame(datasettable[,1:7])


######################################################################
# function - add protein description to blast results from protein annotation at fasta file
######################################################################

#hold_protein_seqs <- 


hold_MAMP_seqs <- data.frame("Protein_Name" = character(0), "MAMP_Hit" = character(0), "Percent_Identity" = character(0),
                             "E-value" = character(0), "MAMP_length" = character(0), "MAMP_Sequence" = character(0), "Genera" = character(0),
                             "Strain_Name" = character(0), "File_Name" = character(0))
  
pb <- txtProgressBar(min = 0, max = length(load_protein_fasta_files), style = 3)

for (i in 1:length(load_protein_fasta_files)){
  
  # read blast resukts and protein fasta file
  read_protein_fasta <-  dss2df(Biostrings::readAAStringSet(load_protein_fasta_files[[i]]))
  read_blast_results <- read.table(file = load_MAMP_blast_result_files[[i]], sep = '\t', header = F)
  colnames(read_blast_results) <- c("Protein_Name","MAMP_Hit","Percent_Identity","E-value","MAMP_length",
                                    "Hit_Start","Hit_End","Hit_Length","Sequence")
  
  # filter proteins in each fasta for those with the same locus tag as the blast results
  grab_right_protein_seq_blast_results <- read_protein_fasta[grepl(paste(read_blast_results$Protein_Name, collapse = "|"), 
                                                                   read_protein_fasta$names),]

  # fix MAMP extension error from blast results
  for (j in 1:nrow(read_blast_results)){
    if (read_blast_results$Hit_Length[j] == read_blast_results$MAMP_length[j]){
      next
    }
    if (read_blast_results$Hit_Length[j] != read_blast_results$MAMP_length[j]){
      subset_protein_seq <- grab_right_protein_seq_blast_results[grepl(read_blast_results$Protein_Name[j], 
                                                                       grab_right_protein_seq_blast_results$names),]
      pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == read_blast_results$MAMP_Hit[j])
      Alignment_between_MAMP_and_Ref <- Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                      subset_protein_seq$seq, type = "global-local", 
                                                                      gapOpening = 100, gapExtension = 100)
      read_blast_results$Sequence[j] <- as.character(Alignment_between_MAMP_and_Ref@subject)
      read_blast_results$Percent_Identity[j] <- Biostrings::pid(Alignment_between_MAMP_and_Ref, type = "PID1")
    }
  }

  # restructure to remove unncessary info
  read_blast_results <- read_blast_results[,c(1,2,3,4,5,9)]
  
  # cross reference with metadata to provide genus, strain, filename info
  get_accession_number <- strsplit(load_protein_fasta_files[[i]], "/")[[1]][6]
  get_accession_number <- paste(strsplit(get_accession_number,"_")[[1]][1],
                                  strsplit(get_accession_number,"_")[[1]][2], sep = "_")
    
  get_strain_info <- subset(datasettable, Assembly_Accession == get_accession_number)
  read_blast_results <- cbind(read_blast_results, 
                              rep(get_strain_info$Genera, nrow(read_blast_results)),
                              rep(get_strain_info$Strain_Name, nrow(read_blast_results)),
                              rep(get_strain_info$Filename, nrow(read_blast_results)))
  
  
  colnames(read_blast_results) <- c("Protein_Name","MAMP_Hit","Percent_Identity","E-value","MAMP_length",
                                    "MAMP_Sequence", "Genera", "Strain_Name", "File_Name")
    
  
  hold_MAMP_seqs <- rbind(hold_MAMP_seqs, read_blast_results)
  setTxtProgressBar(pb, i)
}

close(pb)




##############################################
# Load colors - Set tip labels to match other figures
##############################################

plot_points_Percent_Ident_of_MAMPs <- function(MAMP_of_interest, name_of_MAMP){
  new_plot <- ggplot(MAMP_of_interest, aes(x = Genera, y = Percent_Identity, colour = Genera, fill = Genera)) +
    geom_violin(color = "black") +
    #geom_boxplot(color = "black", outlier.alpha = 0) +
    #geom_quasirandom(varwidth = TRUE, method = "pseudorandom") +
    #geom_jitter(alpha = 0.6, width = 0.2, height = 0.1) +
    scale_color_manual("Genera", values = Genera_colors) +
    scale_fill_manual("Genera", values = Genera_colors) +
    my_ggplot_theme +
    xlab("\n Genera") +
    ylab("Percent Identity \n") +
    scale_y_continuous(limits = c(20,100),
                       breaks = c(20,30,40,50,60,70,80,90,100)) +
    ggtitle(name_of_MAMP) +
    theme(axis.text.x = element_text(angle = 90),
          plot.title = element_text(face = "bold", size =12),
          legend.position = "none") 
  
  return(new_plot)
}


plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'csp22_consensus'), "csp22 and csp22-like")
plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'elf18_consensus'), "elf18 and elf18-like")
plot_points_Percent_Ident_of_MAMPs(subset(hold_MAMP_seqs, MAMP_Hit == 'flg22_consensus'), "flg22 and flg22-like")

##############################################
# Subset to make weblogos
##############################################

# seperate each group into genera
Clav_only <- subset(hold_MAMP_seqs, Genera == "Clavibacter")
Strep_only <- subset(hold_MAMP_seqs, Genera == "Streptomyces")
Leif_only <- subset(hold_MAMP_seqs, Genera == "Leifsonia")
Pseudo_only <- subset(hold_MAMP_seqs, Genera == "Pseudomonas")
Rals_only <- subset(hold_MAMP_seqs, Genera == "Ralstonia")
Curto_only <- subset(hold_MAMP_seqs, Genera == "Curtobacterium")
Xanth_only <- subset(hold_MAMP_seqs, Genera == "Xanthomonas")
Agro_only <- subset(hold_MAMP_seqs, Genera == "Agrobacterium")
Rhodo_only <- subset(hold_MAMP_seqs, Genera == "Rhodococcus")
Rathayi_only <- subset(hold_MAMP_seqs, Genera == "Rathayibacter")



make_me_a_weblogo <- function(file_in){
  logo <- ggseqlogo::ggseqlogo(file_in, seq_type='aa', method = 'bits') +
    theme(panel.grid = element_blank(),
          legend.position = "none", axis.text.x = element_blank(),
          axis.text.y.left = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 12, vjust = 1),
          axis.title.y = element_text(color = "black", size = 12, vjust = 1),
         axis.ticks.x = element_blank())
  
  return(logo)
}


# weblogos ofr csp22 and csp22-like peptides
make_me_a_weblogo(Clav_only[Clav_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Strep_only[Strep_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Leif_only[Leif_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Curto_only[Curto_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Rhodo_only[Rhodo_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Rathayi_only[Rathayi_only$MAMP_Hit %in% 'csp22_consensus',6]) /

  make_me_a_weblogo(Agro_only[Agro_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Pseudo_only[Pseudo_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Rals_only[Rals_only$MAMP_Hit %in% 'csp22_consensus',6]) /
  make_me_a_weblogo(Xanth_only[Xanth_only$MAMP_Hit %in% 'csp22_consensus',6])
  
# weblogos for elf18 and elf18-like peptides
make_me_a_weblogo(Clav_only[Clav_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Strep_only[Strep_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Leif_only[Leif_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Curto_only[Curto_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Rhodo_only[Rhodo_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Rathayi_only[Rathayi_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  
  make_me_a_weblogo(Agro_only[Agro_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Pseudo_only[Pseudo_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Rals_only[Rals_only$MAMP_Hit %in% 'elf18_consensus',6]) /
  make_me_a_weblogo(Xanth_only[Xanth_only$MAMP_Hit %in% 'elf18_consensus',6])
