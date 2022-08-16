MIT_elf18 <- Biostrings::readAAStringSet(filepath = file.choose())
MIT_elf18 <- dss2df(MIT_elf18)

# convert fasta file into 'alignwed' dataframe
MIT_elf18_matrix <- data.frame(MIT_elf18$seq)
MIT_elf18_matrix <- do.call(rbind, strsplit(MIT_elf18_matrix$MIT_elf18.seq, ""))
MIT_elf18_matrix <- as.data.frame(MIT_elf18_matrix)

# based on this paper, we can calculate MIT - https://academic.oup.com/bioinformatics/article/25/9/1125/204722?login=false
# first, we will determine the frequence of each aa for each position
MIT_elf18_AA_freq <- matrix(data = 0, ncol = 20, nrow = 18)
colnames(MIT_elf18_AA_freq) <- Biostrings::AA_STANDARD
MIT_elf18_AA_freq <- as.data.frame(MIT_elf18_AA_freq)
MIT_elf18_AA_freq <- t(MIT_elf18_AA_freq)
colnames(MIT_elf18_AA_freq) <- as.character(seq(1, 18, by=1))


for (i in 1:ncol(MIT_elf18_matrix)){
  calculate_freq_per_position <- t(data.frame(unclass(table(MIT_elf18_matrix[,i]))))
  rownames(calculate_freq_per_position) <- i
  MIT_elf18_AA_freq[match(colnames(calculate_freq_per_position), rownames(MIT_elf18_AA_freq)),i] <- calculate_freq_per_position
}

MIT_elf18_AA_freq_prob <- MIT_elf18_AA_freq
# translate occurance to fequenace numbers
for (i in 1:ncol(MIT_elf18_AA_freq)){
  for (j in 1:nrow(MIT_elf18_AA_freq)){
    MIT_elf18_AA_freq_prob[j,i] <- MIT_elf18_AA_freq[j,i]/sum(MIT_elf18_AA_freq[,i]) 
  }
}


# set the matrix to all 0s
MIT_elf18_matrix_values <- matrix(data = 0, nrow = 18, ncol = 18)
MIT_elf18_matrix_values <- as.data.frame(MIT_elf18_matrix_values)
colnames(MIT_elf18_matrix_values) <- as.character(seq(1, 18, by=1))

# calculate prob. and fill in for each position
for (i in 1:ncol(MIT_elf18_AA_freq_prob)){
  comb_probabilty <- sum(MIT_elf18_AA_freq_prob)
}



$total_mit += $comb_prob*(log( $comb_prob/($pai{$bases[0]} * $pbj{$bases[1]}) )/log(20) );

library(MASS)
rcw_conservation <- function(x) {
  new <- x
  for ( i in 1:nrow(x) ) {
    pos1 <- x[i,1]
    pos2 <- x[i,2]
    
    
    sum1 <- sum(x$Entropy[x$Pos1 == pos1 | x$Pos2 == pos1])
    sum2 <- sum(x$Entropy[x$Pos1 == pos2 | x$Pos2 == pos2])
    #Perform the calculation
    new[i,3] <- x[i,3]/
      (( (sum1+sum2) - (2*x[i,3]) )/
         ( length(x$Entropy[x$Pos1 == pos2 | x$Pos2 == pos2]) + length(x$Entropy[x$Pos1 == pos1 | x$Pos2 == pos1]) - 2 ))
  }
  return(new)
}

z_score <- function(vect_numbers) {
  out <- c()
  for (i in vect_numbers) {
    out <- append(out, (i - mean(vect_numbers))/sd(vect_numbers) )
  }
  return(out)
}

run_analysis <- function(mit, msa = "NA") {
  
  if ( msa != "NA" ) {
    msa$pos <- 1:length(msa[,1])
    mit <- mit[msa$V1[mit$Pos1] > 7 & msa$V1[mit$Pos2] > 7,]
  }
  
  #row column weighting
  mit <- rcw_conservation(mit)
  mit <- na.omit(mit)
  
  #get z-score
  mit$zscore <- z_score(mit$Entropy)
  
  #assign p-vals
  mit$pvalue <- pnorm(mit$zscore, lower.tail = F)
  
  return(mit)
}
