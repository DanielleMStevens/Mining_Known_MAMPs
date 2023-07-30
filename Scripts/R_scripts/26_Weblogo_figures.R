#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 07/06/2020
# Script Purpose: 
# Inputs: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


##############################################
# Subset to make weblogos
##############################################

# seperate each group into genera
Clav_only <- subset(filtered_hold_MAMP_seqs, Genera == "Clavibacter")
Strep_only <- subset(filtered_hold_MAMP_seqs, Genera == "Streptomyces")
Leif_only <- subset(filtered_hold_MAMP_seqs, Genera == "Leifsonia")
Pseudo_only <- subset(filtered_hold_MAMP_seqs, Genera == "Pseudomonas")
Rals_only <- subset(filtered_hold_MAMP_seqs, Genera == "Ralstonia")
Curto_only <- subset(filtered_hold_MAMP_seqs, Genera == "Curtobacterium")
Xanth_only <- subset(filtered_hold_MAMP_seqs, Genera == "Xanthomonas")
Agro_only <- subset(filtered_hold_MAMP_seqs, Genera == "Agrobacterium")
Rhodo_only <- subset(filtered_hold_MAMP_seqs, Genera == "Rhodococcus")
Rathayi_only <- subset(filtered_hold_MAMP_seqs, Genera == "Rathayibacter")




  
  Only_G_pos <- subset(hold_MAMP_seqs, Gram == "+")
  Only_G_neg <- subset(hold_MAMP_seqs, Gram == "-")
  
  # comparisons between Gram-negative and Gram-positives csp22
  make_me_a_weblogo(Only_G_pos[Only_G_pos$MAMP_Hit %in% 'csp22_consensus',6]) /
    make_me_a_weblogo(Only_G_neg[Only_G_neg$MAMP_Hit %in% 'csp22_consensus',6])
  
  # comparisons between Gram-negative and Gram-positives elf18
  make_me_a_weblogo(Only_G_pos[Only_G_pos$MAMP_Hit %in% 'elf18_consensus',6]) /
    make_me_a_weblogo(Only_G_neg[Only_G_neg$MAMP_Hit %in% 'elf18_consensus',6])
  
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
    make_me_a_weblogo(Clav_only[Clav_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Strep_only[Strep_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Leif_only[Leif_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Curto_only[Curto_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Rhodo_only[Rhodo_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Rathayi_only[Rathayi_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    
    make_me_a_weblogo(Agro_only[Agro_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Pseudo_only[Pseudo_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Rals_only[Rals_only$MAMP_Hit %in% 'elf18_consensus',4]) /
    make_me_a_weblogo(Xanth_only[Xanth_only$MAMP_Hit %in% 'elf18_consensus',4])
  
  # weblogos for elf18 and elf18-like peptides
    
    




#-------------- efl18 immunogenic

elf18_immun <- c("AKAKFDRTKPHVNIGTIG","AKAKFERKKPHVNVGTIG","AKAKFERNKLHVNVGTIG",
                 "AKAKFERTKLHVNVGTIG","AKAKFERTKPHVNIGTIG","AKAKFERTKPHVNVGTIG",
                 "AKAKFLREKLHVNVGTIG","AKEKFDRSLPHCNVGTIG","AKEKFDRSLPHVNVGTIG",
                 "AKEKFERNKPHVNVGTIG","AKEKFERSKPHVNVGTIG","AKEKFERTKPHVNVGTIG",
                 "AKERFDRSLPHVNVGTIG","AKGKFERTKPHVNVGTIG","AKSKFERNKPHVNIGTIG",
                 "AKSKFERNKPHVNVGTIG","AKSKFERTKPHVNIGTIG","AKTKFERTKPHVNVGTIG",
                 "ARAKFLREKLHVNVGTIG","PKTAYVRTKPHLNIGTMG","SKEKFERSKPHVNVGTIG",
                 "SKEKFERTKPHVNVGTIG","SKTAYVRTKPHLNIGTMG","GKAKFERTKPHVNIGTIG")


elf18_nonimmuno <- c("PKTAYLRTKPHLNIGTMG","SKKAYVRTKPHLNIGTMG")

make_me_a_weblogo(elf18_immun) + my_ggplot_theme

make_me_a_weblogo(elf18_nonimmuno) + my_ggplot_theme



#----------


csp22_ROS_pos <- c("AQGTVKWFNAEKGFGFIAPEDG","ATGTVKWFNATKGYGFIQPDDG","STGTVKWFNNEKGFGFIAPDDG",
                   "ANGTVKWFNDAKGFGFISPDEG","AQGTVKWFNAEKGYGFIAVDGG","ATGTVKWFNAEKGFGFIAQEGG",
                   "ETGTVKWFNESKGFGFITPDAG","NTGTVKWFNATKGFGFIQPDNG","QSGTVKWFNDEKGFGFITPESG",
                   "QTGTVKWFNDEKGFGFITPQGG","QTGTVKWFNDEKGFGFITPQSG","IKGQVKWFNESKGFGFITPADG",
                   "IKGSVKWFNESKGFGFITPEDG","MNGTVKWFNDAKGFGFITPESG","PNGTVKWFNDAKGFGFISPEDG",
                   "QSGTVKWFNDAKGFGFITPESG","TTGTVKWFNSTKGFGFIQPDNG","DTGTVKWFNTSKGFGFISRDSG",
                   "ETGTVKFFNTDKGFGFIKPDNG","ETGTVKWFNNAKGFGFICPEGG","ETGTVKWFNNAKGFGFICPESG",
                   "MNGIVKWFNDAKGFGFITPESG","MTGTVKWFNNAKGFGFICPAGG","QSGIVKWFNDAKGFGFITPESG",
                   "ISGVVKWFDVAKGFGFIVPDNG","ITGAVKWFDVAKGFGFIVPDNG","ITGVVKWFDVAKGFGFIVPDNG",
                   "MIGLVKWFSPDKGFGFISPTDG","ENGLVKWFNDAKGFGFISRENG","ENGVVKWFNDAKGFGFISRENG",
                   "LNGKVKWFNNAKGYGFIIEDGK","ATGTVKWFNNEKGFGFIAPDDG","ENGTVKWFNDAKGFGFISRENG",
                   "ATGTVKFFAQDKGFGFITPDNG","AHGTLTRWNTDRGFGFITPAQP","ATGTVKWFNAEKGFGFIAPDNG")

csp22_ROS_mid <- c("MTGLVKWFDAGKGFGFITPDNG","ASGKVKWFNNAKGYGFINEEGK","LNGKVKWFNNAKGYGFILEDGK",
                   "PTGTVKWFDHSKGFGFIHPDDG","TQGSVKWFNGEKGFGFIEQDGG","ANGTVKWFNGEKGFGFITVDAV")

#not including cp-46, 47, 63 since not from a annotated csp    
csp22_ROS_neg <- c("PTGKVKWFNSEKGFGFLSRDDG",
                   "YQGRLSDWNDHKGFGFVTPHGG",
                   "YQGRLSDWNDHKGFGFVTPNGG",
                   "MNGTITTWFKDKGFGFIKDENG",
                   "YQGRLRDWNDHKGVGFATPNGG",
                   "FNGIVKNFDLEKGYGFIQPTDG",
                   "PTGKVKFYDDDKGFGFITGDDG",
                   "PTGKVKFYDDEKGFGFISTDDG",
                   "PTGKVKFYDDQKGFGFITGDDG",
                   "PTGKVKFYDEEKGFGFISSDDG",
                   "PTGKVKFYDEEKGFGFISTDDG",
                   "PTGKVKFYDDQKGFGFISGDDG",
                   "PTGKVKFYDEEKGFGFISTDEG",
                   "PTGKVKWYDVDKGFGFLSQEEG",
                   "STGKVIRFDEFKGYGFVAPDEG",
                   "KTGKILRFDEVRGYGFIVPNEG",
                   "PSGRIIKWMTDRGFGFIQEDGA",
                   "YQGRLSDWDDHKGFGFVVPHGG",
                   "QSGEIVDWNDARGFGFIVAAGN",
                   "PTGKVKFYDEDKGFGFISSDDG")



#---------- new version based on clade structure ------------------


csp22_ROS_clade1 <- c("PTGKVKWFNSEKGFGFLSRDDG","PTGKVKFYDDDKGFGFITGDDG",
                      "PTGKVKFYDDEKGFGFISTDDG",
                      "PTGKVKFYDDQKGFGFITGDDG","PTGKVKFYDEEKGFGFISSDDG",
                      "PTGKVKFYDEEKGFGFISTDDG","PTGKVKFYDDQKGFGFISGDDG",
                      "PTGKVKFYDEEKGFGFISTDEG","PTGKVKWYDVDKGFGFLSQEEG",
                      "STGKVIRFDEFKGYGFVAPDEG","KTGKILRFDEVRGYGFIVPNEG",
                      "PTGKVKFYDEDKGFGFISSDDG")

csp22_ROS_clade2 <- c("MNGTITTWFKDKGFGFIKDENG","FNGIVKNFDLEKGYGFIQPTDG")

csp22_ROS_clade3 <- c("YQGRLSDWNDHKGFGFVTPHGG","YQGRLSDWNDHKGFGFVTPNGG",
                      "YQGRLRDWNDHKGVGFATPNGG","YQGRLSDWDDHKGFGFVVPHGG","QSGEIVDWNDARGFGFIVAAGN")


(make_me_a_weblogo(csp22_ROS_pos) + my_ggplot_theme + theme(legend.position = "none")) /
(make_me_a_weblogo(csp22_ROS_mid) + my_ggplot_theme + theme(legend.position = "none")) /  
(make_me_a_weblogo(csp22_ROS_clade1) + my_ggplot_theme + theme(legend.position = "none")) /
(make_me_a_weblogo(csp22_ROS_clade2) + my_ggplot_theme + theme(legend.position = "none")) /
make_me_a_weblogo(csp22_ROS_clade3) + my_ggplot_theme + theme(legend.position = "none")

  
  
  