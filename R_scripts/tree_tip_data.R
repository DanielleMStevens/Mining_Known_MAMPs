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

##############################################
# load data which allows us to change accession name to strain name and add in genus info
##############################################


tip_data <- "./../Mining_for_known_MAMPs_genome_accession_info.xlsx" #choose Mining_for_known_MAMPs_genome_accession_info.xlsx
datasettable <- readxl::read_xlsx(tip_data, col_names = T)
datasettable <- as.data.frame(datasettable[,1:7])
