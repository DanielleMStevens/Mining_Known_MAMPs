## Downloading genomes from NCBI

### 1. Download the genomes
  
  We can use ncbi-genome-download to find what accession we can download for each major genus and then download them on the command line. 
  
  The major genera include the following: 
  |Genera|Clavibacter, Leifsonia, Curtobacterium, Streptomyces, Rathayibacter, Rhodococcus, Agrobacterium, Ralstonia, Xanthomonas, Pseudomonas, Pectobacterium, Dickeya, Erwinia |
  
  We will need to install a package which will allow us to easily download the genomes by accession number from NCBI's refseq.
  
  In order to get accessions for each key genus of bacteria

  ```bash
  conda install ncbi-genome-download
  ncbi-genome-download -s refseq -g Agrobacterium --dry-run bacteria
  ```

  
  I have collected all the accession numbers as well as info about each one into two file stored in the Genome_accession_info directory. I then use the accession name (ex. Erwinia amylovora) to quick filter for accessions that are not either plant/agriculturally related. Once all the information was collected and put into a simple text file, the comman below can be ran:
  
    ```
     ncbi-genome-download --assembly-accessions ./Genome_accession_info/Genome_accessions_to_download.txt -p 6 -r 2 -v --flat-output -F genbank,fasta,protein-fasta bacteria
    ```
    
  where,
    
    ```
     -p 6 : downland 6 genomes at a time in parallel
     -r 2 : retry downloading 2x before moving on
     --flat-out: download all the files in the same place (one directory rather than each isolate having a dedicated directory)
     -v : verbose
     -F 'genbank,fasta,protein-fasta' : download genbank, whole genome fasta, and protein fasta associtaed with the accession number
     ```

Move all the download genomes in directories based on their file type/ending (i.e. genbank files in genbank folder).
 
 ## Setting up database and mining for MAMPs
 
 ### 2. Build the MAMP database
 
 In a text file, save the following MAMP sequences (/MAMP_database/MAMP_elicitor_list.fasta):
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

 ### 3. Run all genomes against balst database


We can then go through each protein fasta file and pull out the peptide from the annotation. In this case, we can use a bash loop to blast each file. 

  ```
  for file in *.faa
  do echo "$file"
  blastp -task blastp-short -xdrop_gap_final 1000 -soft_masking false -query $file -db \ ./../../Mining_Known_MAMPs/MAMP_database/v2/MAMP_blast_db -evalue 1e-4 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue slen qstart qend length qseq" -out $file.txt
  done
  ```

### 4. Processing MAMP database

Using Main_script.R, run though the lines below:

  ```
  ######################################################################
  # set path to data
  ######################################################################
  
    #setwd to where repo was cloned and maintained
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))



  ##############################################
  # Load colors - Set tip labels to match other figures
  ##############################################
  
    #make sure to set path to the same place where the figure 
    source("./package_dependencies.R")
    
    source("./CommonFunctions.R")
    

  ##############################################
  # Load colors - Set tip labels to match other figures
  ##############################################
  
    source("./Figure_colors.R")
  
    source("./Theme_ggplot.R")
  

  ##############################################
  # Run through different sections of each script to process data
  ##############################################

    # load data from blast results
    source("./Loading_raw_data.R")
    
    
    # we then will processes the BLAST results such that they will be organize in a data table
    source("./Process_MAMP_BLAST_results.R")
    
    
    # to be through in our search for microbial MAMPs, we will go back through the annotated genes and pull out
    # the MAMPs that we might have missed in the BLAST search (esspecially for flg22)
    source("./Find_MAMPs_by_annotation.R")
    
  
    # mamps by blast is in "hold_MAMP_seq" and mamps found by annotation are in "All_target_by_annotation"
    source("./combine_blast_and_annotation_results.R")


  ```


 ## Assessing genome diveristy and removing redudnacy/clonality 











### 3. Build Protein Trees of Full Length Sequences and their MAMPs

We now can start building protein trees to understand their evolutionary history in respect to the MAMPs they encode for. We will run MAFFT to build our alignment and IQ-tree of make a maximum likelihood tree from the alignment. In each folder of which the fasta file was saved, the below commands were ran (names changed where needed).

  ```
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

  Run Tree_plotting.R to plot protein and main core gene trees. These trees will be outputed as pdf's at a preset size.

# Pulling out Core genes shared between all strains, species, and genera

In order to build a core gene phyloogeny as well as pull out core genes to calculate Tajima's D, we can use roary to authomate a lot of this. But, unforunately, it doesn't take gbff files (only gff3) and so we need to convert the files before running.

### 1. Prep files for Core Gene analysis

First, we can use a script from [Biocode](https://github.com/jorvis/biocode) which is a set of python scripts for file conversion, among other things. The gbff_to_gff3 script was copied and then the below coommands were ran to complete the install:
  
  ```
  conda deactivate # deactivate conda environment before installing python backage (prevents some weird package conflicts)
  
  apt-get install -y python3 python3-pip zlib1g-dev libblas-dev liblapack-dev libxml2-
  pip3 install biocode
  ```

We can then create a folder to hold the gff3 files

  ```
  for file in *.gbff
  do
        echo "$file"
        python ./../../Mining_Known_MAMPs/python_scripts/convert_genbank_to_gff3.py -i ./../../Whole_Genomes/genbank/$file -o ./../gff3_files/$file.gff3
  done
  ```

  ```
  conda config --add channels r
  conda config --add channels defaults
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda install roary
  ```

# Assessing similarity of genomes in repsect to MAMP diversification

One concern is despite trying to sample for a semi-large data set of genomes in a semi-random manner, if all the genomes are clonal in nature, it's difficult to make any claims in MAMP diveristy. To assess this, we are going to take two different appraoches which should bring us to similar conclusions. One includes building a core gene phylogency (see above) and another is running an ANI analysis. Since I've already built a script to do the same thing for another project, we are going to repurpose that script to do the same in this dataset as well as some downstream analyses. FastANI, software repo seen here, can be used to make these calculations in an automated manner. 


  ```
  ls > MAMP_diversification_ANI_file.txt
  fastANI --rl /path/to/Contig_paths_for_ANI.txt --ql /path/to/Contig_paths_for_ANI.txt -o Clavibacter_ANI_comparison
  ```

 