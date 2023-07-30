


import_pirate_clavibacter <- as.data.table(read_tsv(file = file.choose()))
import_pirate_clavibacter <- as.data.table(read.delim(file = file.choose(), header=F))

colnames(import_pirate_clavibacter) <- c()