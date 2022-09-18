# Table of Contents

- [Packages needed to run this pipeline](#packages-needed-to-run-this-pipeline)
- [Downloading genomes from NCBI](#downloading-genomes-from-ncbi)
  - [1. Download the genomes](#1.-download-the-genomes)
- [Setting up database and mining for MAMPs](#setting-up-database-and-mining-for-mamps)
  - [2. Build the MAMP database](#2.-build-the-mamp-database)
  - [3. Run all genomes against blast database](#3.-run-all-genomes-against-blast-database)
  - [4. Processing Data to Form the MAMP database](#4.-processing-data-to-form-the-mamp-database)
- [Assessing genome diveristy and removing redudnacy/clonality](#assessing-genome-diveristy-and-removing-redudnacy/clonality)
  - [5. Prep genomes for fastANI and run fastANI](#5.-prep-genomes-for-fastani-and-run-fastani)
  - [6. Use fastANI output to filter list and plot dataset diversity](#6.-use-fastani-output-to-filter-list-and-plot-dataset-diversity)


## Packages needed to run this pipeline

Before running all the downstream analyses, we can set up a conda environment with the necessarily software. Below is a brief description of what each package is for.

| Package | Usage | Guide | Citation |
|---------|-------|-------|----------|
| ncbi-genome-download | Will allow us to easily download the genomes by accession number from NCBI's refseq. | [Github Page](https://github.com/kblin/ncbi-genome-download) | N/A |
| fastANI | Calculates whole genome average nucleotide identity at-scale in an all-by-all manner | [Github Page](https://github.com/ParBLiSS/FastANI) | [Paper Link](https://www.nature.com/articles/s41467-018-07641-9) |
| GToTree | Builds phylogenetic trees from whole genomes on the fly based on prepared gene sets | [Github Page](https://github.com/AstrobioMike/GToTree) | [Paper Link](https://academic.oup.com/bioinformatics/article/35/20/4162/5378708) |


```
# Crrates the conda environment
conda create --name myenv

# Uses Conda/Bioconda to install packages
conda install -c bioconda ncbi-genome-download 
conda install -c bioconda fastani 
conda install -c conda-forge -c bioconda -c defaults -c astrobiomike gtotree
conda install -c bioconda blast 

# Activate environment with loaded packages
conda activate myenv
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

  
I have collected all the accession numbers for each genome. Using the accession name (ex. Erwinia amylovora), I quickly filter any accessions that are not either plant/agriculturally related. These accession numbers are listed in a text file in /Analyes/Genome_accession_info directory as well as in a database which stores extra info (Mining_for_known_MAMPs_genome_accession_info.xlsx). The text file list can be used to query and download the accessions using the command below:


```
# Note, the genomes will be downloaded in the local directory of which this command is ran. If the path is changed from the main github repo path, the path to the text file must also be altered.

ncbi-genome-download --assembly-accessions ./Analyses/Genome_accession_info/Genome_accessions_to_download.txt -p 6 -r 2 \\
-v --flat-output -F genbank,fasta,protein-fasta bacteria
```
    
where,
    
```
-p 6 : downland 6 genomes at a time in parallel
-r 2 : retry downloading 2x before moving on
--flat-out: download all the files in the same place (one directory rather 
            than each isolate having a dedicated directory)
-v : verbose
-F 'genbank,fasta,protein-fasta' : download genbank, whole genome fasta, 
    and protein fasta associtaed with the accession number
```

Move all the download genomes in directories based on their file type/ending (i.e. genbank files in genbank folder).

 
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
GSFYSLYFLKDQILNGVNSGHR
```

This fasta file can be used to build a database to use blast to find if anything in the genome shares these sequences. To build the blast database, the below command was ran. Also this should be ran in the same folder as /MAMP_database/MAMP_elicitor_list.fasta. 


```
# to make blast db
makeblastdb -in MAMP_elicitor_list.fasta -parse_seqids -dbtype 'prot' -out MAMP_blast_db
```

### 3. Run all genomes against blast database


We can then go through each protein fasta file and pull out the peptide from the annotation. In this case, we can use a bash loop to blast each file. 

  ```
  for file in *.faa
    do echo "$file"
    blastp -task blastp-short -xdrop_gap_final 1000 -soft_masking false -query $file \\
    -db ./../../Mining_Known_MAMPs/MAMP_database/MAMP_blast_db -evalue 1e-4 -num_threads 12 \\
    -outfmt "6 qseqid sseqid pident evalue slen qstart qend length qseq" -out $file.txt
  done
  ```
  
  where,
  
  ```
  -task blastp-short : specialized version of blastp for short sequences
  -xdrop_gap_final 1000 : heavily penalize dropping mis-hits on edges of MAMP sequences
  -soft_masking false :
  ```


### 4. Processing Data to Form the MAMP database

With the intial serach for MAMPs complete, we will now A) clean up the data by removing any partial hits, B) cross referense the hits by annotation per genome and fill in any missing hits by using local-alignment to the protein with the MAMP of interest to pull out the variant sequenee. 

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

Since we are trying to assess the sequence diversity that has naturally accumulated over time and how that has affected MAMP functionality, we need to do some filtering for clonal isolates so no one sequence is over represented. To do so, we can run all the whole genomes sequences on fastANI to calculate all-by-all ANI values. fastANI tends to inflate values, so we're going to put some strict cutoffs: to be considered clonal, two genomes need to be over 99.99 percent similar and carry the same MAMP eptitope sequnences. Those which are considered clonal by these considerations, duplicates will be removed.

### 5. Prep genomes for fastANI and run fastANI

First we will run an R script designed to pull all the genomes accession number for a particular genus as well as their file path and output them into a directory to run fastANI on each file.

  ```
##############################################
# comparing similarity of genomes to filter for clonality
##############################################

  # write out files to run ANI - run once to generate file to run through fastANI
  source("./10_Parse_genomes_for_ANI_analysis.R")
    
  ```
  
  We can then run the bash command below. Since we have a considerable amount of genomes, this analysis should be ran on a computer that has **no less than 64 Gb of RAM** and at least 10 threads (unless you want to wait many days for it to complete). 
  
  ```bash 
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
  
  
    # ANI Figure  - takes a long time to run so only run once.
    source("./12_ANI_plots.R")
  
  
    # finalizing MAMP list
    source("./13_Finalize_MAMP_list_post_ANI.R")
    ```
  
  
## Phylogentic tree and core gene comparison in respect to MAMPs across species and genera

In order to plot 
  
  ```
  ❯ printf '%s\n' "$PWD"/* >filenames.txt
  #mannually remove files from the path list (based on this file)
  ❯ GToTree -g filenames.txt -H Bacteria -n 4 -j 6 -k -T IQ-TREE
  ```
  
In order to build a core gene phyloogeny as well as pull out core genes to calculate Tajima's D, we can use roary to authomate a lot of this. But, unforunately, it doesn't take gbff files (only gff3) and so we need to convert the files before running.



################################################################################################################ Old code ignore for now




### 6. Prep files for Core Gene analysis

First, we can use a script from [Biocode](https://github.com/jorvis/biocode) which is a set of python scripts for file conversion, among other things. The gbff_to_gff3 script was copied and then the below coommands were ran to complete the install:
  
  ```
  conda deactivate # deactivate conda environment before installing python backage (prevents some weird package conflicts)
  
  apt-get install -y python3 python3-pip zlib1g-dev libblas-dev liblapack-dev libxml2-
  pip3 install biocode
  ```

We can then create a folder to hold the gff files.

  ```
  for file in *.gbff
  do
        echo "$file"
        python ./../../Mining_Known_MAMPs/python_scripts/convert_genbank_to_gff3.py -i ./../../Whole_Genomes_v2/genbank/$file -o ./../gff/$file.gff
  done
  ```
  
But there are subtile variatoins in gff files

  ```
  conda install -c bioconda agat
  for file in *.gff
  do
        echo "$file"
        agat_convert_sp_gxf2gxf.pl -g $file -o ${file}_fixed.gff -gff_version_output 3  
  done


  # move the files into a new folder
  mkdir gff_for_pirate
  mv *.gff_fixed.gff ./../gff_for_pirate
  
  # PIRATE package
  conda install pirate 

  # optional dependencies for plotting figures in R
  conda install r==3.5.1 r-ggplot2==3.1.0 r-dplyr==0.7.6 bioconductor-ggtree==1.14.4 r-phangorn==2.4.0 r-gridextra
  conda install -c bioconda diamond
  
  # run pirate locally in new gff folder
  PIRATE -i ~/Documents/Mining_MAMPs/Whole_Genomes_v2/gff_for_pirate/ -s "70,90" -k "--cd-step 2 --cd-low 90" -a -r -t 14 
  ```
  
  We then use datasettable and genomes_to_check in the R-console to view the dataframe and remove gff files manually which are considered clonal by previous analyses. Once we 

  ```
  # run the below with the name of the environment, run once
  #conda create --name roary_analysis

  
  
  


### 3. Build Protein Trees of Full Length Sequences and their MAMPs

We now can start building protein trees to understand their evolutionary history in respect to the MAMPs they encode for. We will run MAFFT to build our alignment and IQ-tree of make a maximum likelihood tree from the alignment. In each folder of which the fasta file was saved, the below commands were ran (names changed where needed).

  ```
  # for MAMP sequnece trees:
  ❯ mafft --reorder --thread 12 --maxiterate 1000 csp22.fasta > "csp22_alignment"
  
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
  
  
  # MAMP trees 
  ❯ mafft --reorder --thread 12 --maxiterate 1000 --localpair --op 3 csp22_for_tree.fasta > "csp22_MAMP_alignment"
  iqtree -s csp22_MAMP_alignemnt -bb 1000 -T AUTO -v -m TEST
  ```

  Run Tree_plotting.R to plot protein and main core gene trees. These trees will be outputed as pdf's at a preset size.

# Assessing similarity of genomes in repsect to MAMP diversification

One concern is despite trying to sample for a semi-large data set of genomes in a semi-random manner, if all the genomes are clonal in nature, it's difficult to make any claims in MAMP diveristy. To assess this, we are going to take two different appraoches which should bring us to similar conclusions. One includes building a core gene phylogency (see above) and another is running an ANI analysis. Since I've already built a script to do the same thing for another project, we are going to repurpose that script to do the same in this dataset as well as some downstream analyses. FastANI, software repo seen here, can be used to make these calculations in an automated manner. 


  ```
  ls > MAMP_diversification_ANI_file.txt
  fastANI --rl /path/to/Contig_paths_for_ANI.txt --ql /path/to/Contig_paths_for_ANI.txt -o Clavibacter_ANI_comparison
  ```




 
