hold_cspD <- c("MQNGKVKWFNNEKGFGFIEVEGGDDVFVHFTAIEGDGYKSLEEGQEVSFEIVEGNRGPQASNVVKL")
hold_cspB <- c("MLEGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTKEA")

test2 <- Biostrings::pairwiseAlignment(hold_cspD, All_target_by_annotation$seq[4], type = "global")
test2
Biostrings::pid(test2, type = "PID1")
