# Table of Contents

- [Packages needed to run this pipeline](#packages-needed-to-run-this-pipeline)

- [Downloading genomes from NCBI](#downloading-genomes-from-ncbi)
    - [1. Download the genomes](#1-download-the-genomes)

- [Setting up database and mining for MAMPs](#setting-up-database-and-mining-for-mamps)
    - [2. Build the MAMP database](#2-build-the-mamp-database)
    - [3. Run all genomes against blast database](#3-run-all-genomes-against-blast-database)
    - [4. Processing Data to Form the MAMP database](#4-processing-data-to-form-the-mamp-database)
  
- [Assessing genome diveristy and removing redudnacy/clonality](#assessing-genome-diveristy-and-removing-redudnacy/clonality)
    - [5. Prep genomes for fastANI and run fastANI](#5-prep-genomes-for-fastani-and-run-fastani)
    - [6. Use fastANI output to filter list and plot dataset diversity](#6-use-fastani-output-to-filter-list-and-plot-dataset-diversity)
    - [7. Phylogentic tree in respect to MAMPs abundance across genera](#7-phylogentic-tree-in-respect-to-mamps-abundance-across-genera)

- [Characterizing Diversity of MAMPs and their MAMP-Encoded Proteins](#characterizing-diversity-of-mamps-and-their-mamp-encoded-proteins)
    - [8. Displaying the relatedness of epitope variants in respect to immmunogenicity](#8-displaying-the-relatedness-of-epitope-variants-in-respect-to-immmunogenicity)
    
- [Evolution of the Gene MAMPs are Encoded on](#evolution-of-the-gene-mamps-are-encoded-on)
    - [9. Assessing immunogenicity outcomes following diversifcation of a common ancestor](#9-assessing-immunogenicity-outcomes-following-diversifcation-of-a-common-ancestor)
    - [10. Evolution of CSPs](10-evolution-of-csps)
    - [11. Core gene analyses to assess selection (Tajma's D)](#11-core-gene-analyses-to-assess-selection-(tajma's-d))


<br>


## Packages needed to run this pipeline

Before running all the downstream analyses, we can set up a conda environment with the necessarily software. Below is a brief description of what each package is for.

| Package | Usage | Guide | Citation |
|---------|-------|-------|----------|
| ncbi-genome-download | Will allow us to easily download the genomes by accession number from NCBI's refseq. | [Github Page](https://github.com/kblin/ncbi-genome-download) | N/A |
| fastANI | Calculates whole genome average nucleotide identity at-scale in an all-by-all manner | [Github Page](https://github.com/ParBLiSS/FastANI) | [Paper Link](https://www.nature.com/articles/s41467-018-07641-9) |
| GToTree | Builds phylogenetic trees from whole genomes on the fly based on prepared gene sets | [Github Page](https://github.com/AstrobioMike/GToTree) | [Paper Link](https://academic.oup.com/bioinformatics/article/35/20/4162/5378708) |
| Blast+ | Enables running blast on the command line | [NCBI Page](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) | N/A |
| Agat | A tool which helps convert between various iterations of Gff file formats | [Github Link](https://github.com/NBISweden/AGAT) | N/A - Zenodo DOI - see Github Page |
| Pirate | Determines Core genes based on a group of genomes | [Github Page](https://github.com/SionBayliss/PIRATE) | [Paper Link](https://academic.oup.com/gigascience/article/8/10/giz119/5584409) |
| FastTree | FastTree infers approximately-maximum-likelihood phylogenetic trees | [Github Page](tbd) | [Paper Link](TBD) |
| HMMER | Finds sequence homologs using profile hidden Markov models | [Github Page](https://github.com/EddyRivasLab/hmmer) | N/A |
| MMseqs2 | ultra fast and sensitive sequence search and clustering suite | [Github Page](https://github.com/soedinglab/MMseqs2) | [Paper Link](https://academic.oup.com/bioinformatics/article/37/18/3029/6178277?login=true) |
| TrimAL | Tool to trim alignments | [Github Page](https://github.com/inab/trimal) | [Paper Link](http://trimal.cgenomics.org/_media/trimal.2009.pdf) |


<br>

```
# Creates the conda environment 
conda create --name mining_mamps

# Uses Conda/Bioconda to install packages
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -c bioconda ncbi-genome-download 
conda install -c bioconda fastani 
conda install -c conda-forge -c bioconda -c defaults -c astrobiomike gtotree
conda install -c bioconda blast 
conda install -c bioconda agat
conda install -c bioconda pirate 
conda install r==3.5.1 r-ggplot2==3.1.0 r-dplyr==0.7.6 bioconductor-ggtree==1.14.4 r-phangorn==2.4.0 r-gridextra
conda install -c bioconda diamond
conda install -c bioconda fasttree
conda install -c bioconda hmmer
conda install -c conda-forge -c bioconda mmseqs2
conda install -c bioconda mafft
conda install -c bioconda iqtree
conda install -c bioconda meme


# Activate environment with loaded packages
conda activate mining_mamps
```

## Downloading genomes from NCBI

### 1. Download the genomes

We can use ncbi-genome-download to find which genome accessions we can download for each genus and then download them on the command line. 

We decided to focus on several genera that have many known phytopathogens as well as all genera from actinobacterial pathogens. One aspect to note is that not all genomes are necessarily from pathogens as **(1)** this would be difficult to verify at scale and **(2)** previous work has shown examples of symbionts interacting similarily as pathogens in the content of MAMPs and PRR interactions. Therefore, any genome which is plant-, agriculture-, (and in some cases soil-) associated is included in the dataset.

| Gram-type | Genera|
| ------------- | --------------------|  
| Gram-positive | Clavibacter, Leifsonia, Curtobacterium, Streptomyces, Rathayibacter, Rhodococcus |
| Gram-negative | Agrobacterium, Ralstonia, Xanthomonas, Pseudomonas, Pectobacterium, Dickeya, Erwinia |

  
In order to get accessions for each key genus of bacteria

```
ncbi-genome-download -s refseq -g Agrobacterium --dry-run bacteria
```

  
I have collected all the accession numbers for each genome. Using the accession name (ex. Erwinia amylovora), I filter any accessions that are not plant/agriculturally related. These accession numbers are listed in a text file in /Analyses/Genome_accession_info directory as well as in a database which stores extra info (Mining_for_known_MAMPs_genome_accession_info.xlsx). The text file list can be used to query and download the accessions using the command below:


```
# Note, the genomes will be downloaded in the local directory of which this command is ran. 
# If the path is changed from the main github repo path, the path to the text file must also be altered.

ncbi-genome-download --assembly-accessions ./Analyses/Genome_accession_info/Genome_accessions_to_download.txt -p 6 \\
-r 2 -v --flat-output -F genbank,fasta,protein-fasta bacteria
```
    
where,
    
```
-p 6 : downland 6 genomes at a time in parallel
-r 2 : retry downloading 2x before moving on
--flat-out: download all the files in the same place (one directory rather 
            than each isolate having a dedicated directory)
-v : verbose
-F 'genbank,fasta,protein-fasta' : download genbank, whole genome fasta, 
    and protein fasta associated with the accession number
```

I moved all the download genomes in directories based on their file type/ending (i.e. genbank files in genbank folder).

 
## Setting up database and mining for MAMPs
 
### 2. Build the MAMP database
 
In a text file, save the following MAMP sequences (/Analyses/MAMP_blast/MAMP_database/MAMP_elicitor_list.fasta):
 
```
>csp22_consensus
AVGTVKWFNAEKGFGFITPDDG
>elf18_consensus
SKEKFERTKPHVNVGTIG
>flg22_consensus
QRLSTGSRINSAKDDAAGLQIA
>flgII-28
ESTNILQRMRELAVQSRNDSNSATDREA
>nlp20_consensus
AIMYSWYFPKDSPVTGLGHR
```

This fasta file can be used to build a database to use blast to find if anything in the genome shares these sequences. To build the blast database, the below command was ran. Also this should be ran in the same folder as /MAMP_database/MAMP_elicitor_list.fasta. 


```
# to make blast db
makeblastdb -in MAMP_elicitor_list.fasta -parse_seqids -dbtype 'prot' -out MAMP_blast_db
```

### 3. Run all genomes against blast database


We can then go through each protein fasta file and pull out the possible peptide. In this case, we can use a bash loop to blast each file. 

  ```
  for file in *.faa
    do echo "$file"
    blastp -task blastp-short -xdrop_gap_final 1000 -soft_masking false -query $file \\
    -db ./../../Mining_Known_MAMPs/Analyses/MAMP_blast/MAMP_database/MAMP_blast_db \\
    -evalue 1e-4 -num_threads 12 \\
    -outfmt "6 qseqid sseqid pident evalue slen qstart qend length qseq" -out $file.txt
  done
  ```
  
  blastp -task blastp-short -xdrop_gap_final 1000 -soft_masking false -query $file -db ./../../Mining_Known_MAMPs/Analyses/MAMP_blast/MAMP_database/MAMP_blast_db -evalue 1e-4 -num_threads 12 -outfmt "6 qseqid sseqid pident evalue sl
en qstart qend length qseq" -out $file.txt
  
  where,
  
  ```
  -task blastp-short : specialized version of blastp for short sequences
  -xdrop_gap_final 1000 : heavily penalize dropping mis-hits on edges of MAMP sequences
  -soft_masking false : to prevent false-negative annotations of sequences of interest
  ```


### 4. Processing Data to Form the MAMP database

With the intial serach for MAMPs complete, we will now A) clean up the data by removing any partial hits, B) cross referense the hits by annotation per genome and fill in any missing hits by using local-alignment to the protein with the MAMP of interest. This will allow us to pull out the variant peptide sequence. 

Using Main_script.R, run though the lines below:

  ```

  ######################################################################
  # set path to data
  ######################################################################

  #setwd to where repo was cloned and maintained
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))


  ##############################################
  # Load packages and functions
  ##############################################

  #make sure to set path to the same place where the figure 
  source("./01_Package_dependencies.R")
  
  source("./02_CommonFunctions.R")
  

  ##############################################
  # Load colors and ggplot theme
  ##############################################
  
  source("./03_Figure_colors.R")
  
  source("./04_Theme_ggplot.R")
  

  ##############################################
  # Load data and Run through different sections of each script to process data MAMP data 
  ##############################################

  # load data from blast results
  source("./05_Loading_raw_data.R")
  
  
  # we then will processes the BLAST results such that they will be organize in a data table
  source("./06_Process_MAMP_BLAST_results.R")
  
  
  # to be through in our search for microbial MAMPs, we will go back through the annotated genes and pull out the 
  # MAMPs that we might have missed in the BLAST search (esspecially for flg22)
  source("./07_Find_MAMPs_by_annotation.R")
  

  # mamps by blast is in "hold_MAMP_seq" and mamps found by annotation are in "All_target_by_annotation"
  source("./08_Combine_blast_and_annotation_results.R")


  # filter for partial proteins and mannual remove a small number of off-targets
  source("./09_Filter_for_partial_proteins.R")

  ```


 ## Assessing genome diveristy and removing redudnacy/clonality 

Since we are trying to assess the sequence diversity that has naturally accumulated over time and how that has affected MAMP functionality, we need to do some filtering for clonal isolates so no one sequence is over represented. To do so, we can run all the whole genomes sequences on fastANI to calculate all-by-all ANI values. fastANI tends to inflate values, so we're going to put some strict cutoffs: to be considered clonal, two genomes need to be over 99.999 percent similar and carry the same MAMP eptitope sequnences. Those which are considered clonal by these considerations, duplicates will be removed. This number was also based on using the Clavibacter genus as a reference since I'm more personally familair with those genomes and if this number was reasonable.

### 5. Prep genomes for fastANI and run fastANI

First we will run an R script designed to pull all the genomes accession number for a particular genus as well as their file path and output them into a directory to run fastANI on each file.

  ```R
##############################################
# comparing similarity of genomes to filter for clonality
##############################################

  # write out files to run ANI - run ONCE to generate file to run through fastANI
  source("./10_Parse_genomes_for_ANI_analysis.R")
    
  ```
  
  We can then run the bash command below. Since we have a considerable amount of genomes, this analysis should be ran on a computer that has **no less than 64 Gb of RAM** and at least 10 threads (unless you want to wait many days for it to complete). 
  
  ``` bash
  for file in *.txt
    do echo "$file"
    fastANI --rl ./$file --ql ./$file -t 14 -o ${file}_ANI_comparison
  done
  ```
  
  For fastANI, the following commands are:
  
  ```
  --rl: path to the file which contains path inforamtion for the genomes
  --ql: the same file and file path for --rl, this allows us to run all-by-all comparisons
  -t 14: 14 threads to run calculations in parallel, adjust this as needed 
  -o: output file name
  ```
  
### 6. Use fastANI output to filter list and plot dataset diversity
  
  The output files can be parsed for their ANI values, compared in respect to the MAMP eptitopes found, and used to filter and create a final MAMP list. Additionally, we will create a figure showing the ANI values as a representative of a diverse dataset.
  
  Using Main_script.R, run though the lines below:
  
    ```
    # parsing ANI values and comparing them to MAMP sequences to remove clonal genomes
    source("./11_ANI_analysis.R")  
    
    
    # finalizing MAMP list
    source("./12_Finalize_MAMP_list_post_ANI.R")
    
    
    # ANI Figure  - takes a long time to run so only run once.
    # Part of Supplemental Figure 1
    source("./13_ANI_plots.R")
    ```
  
  
# Assessing similarity of genomes in repsect to MAMP diversification

One concern is despite trying to sample for a semi-large data set of genomes in a semi-random manner, if all the genomes are clonal in nature, it's difficult to make any claims in MAMP diveristy. To assess this, we are going to take two different appraoches which should bring us to similar conclusions. One includes building a core gene phylogency (see above) and another is running an ANI analysis. Since I've already built a script to do the same thing for another project, we are going to repurpose that script to do the same in this dataset as well as some downstream analyses. FastANI, software repo seen here, can be used to make these calculations in an automated manner. 


  ```
  ls > MAMP_diversification_ANI_file.txt
  fastANI --rl /path/to/Contig_paths_for_ANI.txt --ql /path/to/Contig_paths_for_ANI.txt -o Clavibacter_ANI_comparison
  ```

## Assessing MAMP Abundance across diverse bacteria 

Now that we have a finalized MAMP list across a large number of genomes from diverse agriculturally-associated bacteria, I wanted to see how abundant each MAMP is across the genera assessed. To do so, I want to do a phylogenetic approach which is rhobust but not computationally intensive. To do so, I choose to use GToTree (into here)[https://github.com/AstrobioMike/GToTree], which assess 74 housekeeping genes found in bacteria. This should be more rhobust than a MSLA tree and less intensive than a core gene phylogeny.
  
### 7. Phylogentic tree in respect to MAMPs abundance across genera

In order to plot the distribution of these MAMPs, I first needed to list all the genomes assessed in this study using the command below.
  
  ```
  # First we will print all the gbff genome files into a text file
  printf '%s\n' "$PWD"/* > filenames.txt
  ```
  
  I will then mannually remove genomes listed from the file if there were found as clonal based on fastANI comparisons. This final list in file named filenames.txt will be passed to GToTree with the following below parameters.
  
  ```
  GToTree -g filenames.txt -H Bacteria -n 8 -j 10 -k -G 0.2 -o GToTree_output
  ```
  
  where,
  
  ```
  -H Bacteria : HMM models for bacteria to bulld key genes representative of a bacterial dataset
  -n 8 : Number of cpus to use during hmm search
  -j 10 : Number of jobs to run during parallelizable steps
  -k : Individual protein alignment files will be retained
  -G 0.2 : Genome minimum gene-copy threshold
  ```
  
  Note -X parameter, which does not speed up alignments, was not used for this analyses becuase 1. it seems to be broken on my installation an 2. 100% percision isn't needed for this kind of anaylses. Make sure the tree files are moved to the appropriate directory (Analyses/All_bacteria_phylogenetic_tree/Figure1B_Phylogenetic_tree/). We then will return to the Main_script.R to process the tree data as well as make some plots on general comparisons of the diversity of MAMP eptiopes.
  
    Using Main_script.R, run though the lines below:
  
    ```
    ##############################################
    # Dyanamic changes of MAMPs from diverse bacteria - Figure 1
    ##############################################
      
      #*******************************************
      # go back to methods.md to run phylogeneitc analysis
      #*******************************************
      
      # phylogenetic tree of diverse bacteria with MAMP eptitope number plotted on
      source("./14_Core_tree_with_MAMP_number_plotted.R")
    
      
      # ggplot figures such as violin plots
      # includes Figure 1C, D, E; as well as part of Supplemental Figure 1
      source("./15_Comparison_between_MAMP_and_consensus_ggplot_figures.R")
      
      
      # Assessing variation within each MAMP epitope and in a positional manner 
      # includes Figure 1G
      source("./16_Assessing_MAMP_variation.R")
      
      
      # how many MAMPs are there, the variation that extists, and how many/offen plot
      # includes Supplemental Figure 4
      source("./17_Abundance_of_MAMPs.R")

    ```
  
  We see that MAMPs how different trajectories in how they have changed over time. This is likely due to a combination of factors including differences in selection due to protein function/necessity to life, protein abundance, and cognate receptor conservation. But there are futher questions on how this diversity arrise and what this means in the context of impacting immune perception.
  
  
  
## Characterizing Diversity of MAMPs and their MAMP-Encoded Proteins
  
First, we will want to plot the immune outcomes for each epitope tested (via ROS burst) and then secondarily, we will want to assess if this immune output is driven by the associated MAMP-dervived gene. To assess this, we will 1) make a cladogram based on the sequences of the epitopes tested and 2) make a phylogenic tree and plot information in respect to the genera the protein is derived from, the amino acid similarity, and the immnogencity to understand which and how much different eptiopes are changing and if that is correlated with immunogencity, 

### 8. Displaying the relatedness of epitope variants in respect to immmunogenicity


After completeing the screen, all the epitope sequences for each MAMP were copied into a fasta file (in Analyses/ROS Screen/epitope_tree_for_ROS). We then built the tree using both mafft and FastTREE. For each epitope, run the below code (though altering the file names):


   ```
   # build tree for elf18 variants
   mafft --auto elf18_variants.fasta > elf18_variants_aligned
   FastTree elf18_variants_aligned > elf18_variants_tree
   
   # build tree for csp22 variants
   mafft --auto csp22_variants.fasta > csp22_variants_aligned
   FastTree csp22_variants_aligned > csp22_variants_tree
   ```

The outputs for both elf18 and csp22 variants were saved in the Analyses/ROS Screen/epitope_tree_for_ROS directory. We will then return to the main page and run the following script.

    Using Main_script.R, run though the lines below:
    
    ```
    ##############################################
    # Diverse epitopes and their impact on plant immmune perception - Figure 2 and 3
    ##############################################

    # Creating phylogenic tree with immunogenicity data (via ROS assays)
    source("./18_ROS_Screen_tree.R")
    
    ```

These output trees are used to build Figure 2A and Figure 3A. The tips are colored manually as well as adding the of the data using Inkscape. This script will also import the immunogenicity conclusions which will be mapped onto the phylogenetic trees in #9: Evolution of CSP variants (Analyses/ROS_Screen/ROS_Screen_data.xlsx).

## Evolution of the Gene MAMPs are Encoded on

### 9. Assessing immunogenicity outcomes following diversifcation of a common ancestor

Looking at the cladograms, there seems to be a trend of which epitopes are related and how that impacts immune perception. But keep in mind, the tree were only built off the epitope sequence not the whole protein. To assesss if the MAMP-derived variants which have different immunological outcomes have diverved from a common ancestor or potentially via convergent evolution, we will built a whole protein evolutionary tree (rather than just the epitope). With EF-Tu being so conserved as a housekeeping gene, this will relatively trival. However, with CSPs being so different in sequence, yet conserved in structure, this will be a little trickier. For the latter, we will try a similar approach to what Ksenia's lab has done with NLR work ([link](https://github.com/krasileva-group/hvNLR)) where we will build the tree based on a conserved domain and plot additional informaiton onto the tree.

First, we will need to re-pull out all the protein sequences for each MAMP variant. We can do this using the main script: 


    Using Main_script.R, run though the lines below:
    
    ```
    ##############################################
    # Diverse epitopes and their impact on plant immmune perception - Figure 2 and 3
    ##############################################

    # Re-pull whole protein sequences for MAMP hits and write to fasta file
    source("./19_Write_MAMP_hits_to_fasta.R")
  
    ```
    
For EF-Tu, we moved the fasta file to the following directory (Analyses/Protein_alignments_and_trees/EfTu) and then run the below code:

  ``` 
  # build tree for EF-Tu homologs
  
  # for full length proteins, 
  mafft --reorder --thread 12 --maxiterate 1000 --localpair EFTu_full_length.fasta > "EFTu_full_length_alignment"
  
  # --localpair, slowest but most accurate method of alignment
  # --reorder, reorder entries in fasta file to improve alignment
  
  
  iqtree -s EFTu_full_length_alignment -m MFP -bb 1000 -T AUTO -v
  # -s, input alignment file
  # -st, file type (in this case amino acids, hence AA)
  # -bb 1000, number of ultrafast bootstrapping ran on the tree
  # -mtree, iterate thorugh all models to find the best one
  # -nt 12, number of threads used to run analysis (I have max 16)

  ```


For CSPs, we first need to obtain the conserved CSP domain model and use HMMER to extract this domain from each CSP sequence. This was previously found on the Pfam website but as of January 2023, that website has been deprcipated. Details of the domain can be now found on the InterPro site ([link here](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR011129/)). A copy of this model was downloaded and move into a directory in this Github repo (Analyses/Characterizing_CSPs/Hmm_models/CSD_model).


  ```
  # Search all the protein fasta files for the CSP doamin
  hmmsearch -A CSD_alignment.stk --tblout csp_domains.txt -E 1 --domE 1 --incE 0.01 --incdomE 0.04 --cpu 8 \//
  ./Hmm_modles/CSD.hmm csp_full_length.fasta 
  
  # Convert the output from hmmersearch into a fasta file
  esl-reformat fasta CSD_alignment.stk > reformat_CSD_hits.fasta
  ```

  However, before we can build a tree off of this, we need to remove redundancies and clean up header names. Based off of previous manual inspecion
  
  
  We will then build our alignment and phylogenetic tree using the following commands.
  
    ``` 
    # build tree for CSPs homologs and paralogs
  
    # for full length proteins, 
    mafft --reorder --thread 12 --maxiterate 1000 --localpair reformat_CSD_hits.fasta > "full_tree_csp_domain_aligned"
  
    # --localpair, slowest but most accurate method of alignment
    # --reorder, reorder entries in fasta file to improve alignment
  
    iqtree -s full_tree_csp_domain_aligned -st AA -bb 1000 -T AUTO -m MFP -keep-ident-safe
  
    
    # -s, input alignment file
    # -st, file type (in this case amino acids, hence AA)
    # -bb 1000, number of ultrafast bootstrapping ran on the tree
    # -mtree, iterate thorugh all models to find the best one
    # -nt 12, number of threads used to run analysis (I have max 16)
    ```
  
This tree may take awhlie o build as it is quite large. Once it is completed, we will build our phylogenetic trees using genera, MAMP epitopes, and immunogenicity outcome data to better understand the evolution of these proteins.  

```
  # for full length proteins, 
  mafft --reorder --thread 12 --maxiterate 1000 --localpair csp22_full_length.fasta > "csp22_full_length_alignment"
  # --localpair, slowest but most accurate method of alignment
  # --reorder, reorder entries in fasta file to improve alignment
  
  iqtree -s csp22_full_length_alignment -st AA -bb 1000 -mtree -nt 12 -keep-ident	-safe
  # -s, input alignment file
  # -st, file type (in this case amino acids, hence AA)
  # -bb 1000, number of ultrafast bootstrapping ran on the tree
  # -mtree, iterate thorugh all models to find the best one
  # -nt 12, number of threads used to run analysis (I have max 16)
```

### 10. Evolution of CSPs

Considering each MAMP has a different evolutionary trajectory, we then wanted to better understand the diversity and evolution of the genes/proteins of which the MAMPs are encoded. This is of particular interest for CSPs as they are so diverse and so many copies are present. To do so, we will use a variety of techniques including the phylogenetic tree built in #8 as well as mmseqs2 and MEME suite. First, we will use mmseq2 to cluster



Details for catagotizing 

    ```
    # Re-pull whole protein sequences for MAMP hits and write to fasta file
    source("./23_by_genera_csps_to_fasta.R")
    
    # for each outputed genera specific fasta file, we will first catagorize using mmseq 2 and the 
    # outputs will go into directories labeled by each genera
    ❯ mmseqs easy-cluster Rathayibacter_CSPs.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1
    
    ```

We will then build phylogenetic trees from the csp domain. To do this, we need to extra the csp domain, reformat the flie, align the sequences and build the tree.


    ```
    # For each genera sepecifc csp fasta file, run the following commands.
    hmmsearch -A CSD_alignment.stk --tblout csp_domains.txt -E 1 --domE 1 --incE 0.01 --incdomE 0.04 --cpu 8 ./Hmm_modles/CSD.hmm csp_full_length.fasta 
  
    # Convert the output from hmmersearch into a fasta file
    esl-reformat fasta CSD_alignment.stk > reformat_CSD_hits.fasta
    
    ```
    
When we combine these data together, it reveals a clear strategy of which CSPs are related to each other. We can then cross-reference these groups by using all-by-all blastp search.


## Determining Core genes to assess selection of MAMP-endcoded genes compared to other conserved genes

In MAMP-focused resarch, there is an idea that MAMPs undergo co-evolution with their cognant PRR receptors. Therefore, in order to evade perception, MAMPs are thought to be under positive selection compared to other gene counterparts (which are typically under negative selection). In order to assess this on an epitope and overall gene basis, we need to first determine which genes are considered 'core' accross all the genomes. However, determining which genes are considered 'core' accross genera so diverse may be difficult and computationally rigorus. Therefore, we aim to find intra-genus core genes and calculate a variety of pop. gene stats.
  

### 11. Core gene analyses to assess selection (Tajma's D)

We will use Pirate to determine which genes are considered 'core' for each genus of bacteria but, unforunately, it doesn't take gbff files (only gff3) and so we need to convert the files before running. First, we moved the list of accessions into new text files. We will then remove those that are redudant based on previous ANI analyses as similar to before. Then we will temperarily re-download the gbff files for each genus into dedicated folders (per genus) and convert them into gff3 files (thus removing the gbff files).

For each text file in each genus specific folder:

  ```
  ncbi-genome-download --assembly-accessions ./Clavibacter_genome_list.txt -p 6 -r 2 -v --flat-output -F genbank bacteria
  cds-fasta
  ```
  
For each genus, the gbff files were unzip and convert into gff files. But even then those files needed
#https://pypi.org/project/biocode/

  ```
  for file in *.gbff
  do
        echo "$file"
        python ./../../Mining_Known_MAMPs/python_scripts/convert_genbank_to_gff3.py -i ./../../Whole_Genomes_v2/genbank/$file -o ./../gff/$file.gff
  done
  ```







But there are subtile variatoins in gff files

  ```
  for file in *.gff
  do
        echo "$file"
        agat_convert_sp_gxf2gxf.pl -g $file -o ${file}_fixed.gff -gff_version_output 3  
  done


  # move the files into a new folder
  mkdir gff_for_pirate
  mv *.gff_fixed.gff ./../gff_for_pirate
  

  
  # run pirate locally in new gff folder
  PIRATE -i ~/Documents/Mining_MAMPs/Whole_Genomes_v2/gff_for_pirate/ -s "70,90" -k "--cd-step 2 --cd-low 90" -a -r -t 14 
  ```
  


  
  
  


################################################################################################################ Old code ignore for now




### 6. Prep files for Core Gene analysis

First, we can use a script from [Biocode](https://github.com/jorvis/biocode) which is a set of python scripts for file conversion, among other things. The gbff_to_gff3 script was copied and then the below coommands were ran to complete the install:
  
  ```
  conda deactivate # deactivate conda environment before installing python backage (prevents some weird package conflicts)
  
  apt-get install -y python3 python3-pip zlib1g-dev libblas-dev liblapack-dev libxml2-
  pip3 install biocode
  ```

We can then create a folder to hold the gff files.


  
  We then use datasettable and genomes_to_check in the R-console to view the dataframe and remove gff files manually which are considered clonal by previous analyses. Once we 

  ```
  # run the below with the name of the environment, run once
  #conda create --name roary_analysis
  ```
  
  
  


 
