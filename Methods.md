## Downloading genomes from NCBI

### 1. Download Conda

  We can use conda to help download software needed to find the orthologs we desire and processes them for phylogenetic analsyes.

    # go to this website https://www.anaconda.com/distribution/
    # select the 3rd distribution if you don't already use anaconda 
    bash Anaconda-latest-Linux-x86_64.sh


    # run so your computer can find conda with your path
    #for anaconda 2 :

    export PATH=~/anaconda2/bin:$PATH

    #for anaconda 3 :

    export PATH=~/anaconda3/bin:$PATH


  if you need to initalize conda, write the command:
    
    $ conda init {shell script - aka bash}
  
  then restart ternminal, it should look like this: 
    
    $ (base) danimstevens@danimstevens-MS-7B93:

  to set up other conda channels
  
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

  After conda is initalized, we can use it to stabaily download a variety of bioinformatic software. Search this website to see what they have available: https://anaconda.org/bioconda


### 2. Download the genomes
  
  We can use ncbi-genome-download to find what accession we can download for each major genus and then download them on the command line. The major 6 genera include Clavibacter, Leifsonia, Curtobacterium, Streptomyces, Rathayibacter, and Rhodococcus. I have collected all the accession numbers as well as info about each one into two file stored in the Genome_accession_info directory. First we will need to install a package which will allow us to easily download the genomes by accession number from NCBI's refseq.
  
    conda install ncbi-genome-download
    
  Once that has installed, run the following command to download the genomes we are interested in. 

In order to get accessions for each key genus of bacteria

```bash
ncbi-genome-download -s refseq -g Agrobacterium --dry-run bacteria
```

Once all the information is collected and put into a simple text file (/Genome_accession_info/Genome_accessions_to_download.txt), the comman below can be ran:


```bash
ncbi-genome-download --assembly-accessions ./Genome_accession_info/Genome_accessions_to_download.txt \ 
-p 4 -r 2 -v --flat-output -F genbank,fasta,protein-fasta bacteria
```

where,

```
-p 4 : downland 4 genomes at a time in parallel
 -r 2 : retry downloading 2x before moving on
 --flat-out: download all the files in the same place (one directory rather than each isolate having a dedicated directory)
 -v : verbose
 -F 'genbank,fasta,protein-fasta' : download genbank, whole genome fasta, and protein fasta associtaed with the accession number
 ```
 Move all the download genomes in directories based on their file type/ending (i.e. genbank files in genbank folder).
 
 ## Setting up database and mining for MAMPs
 
 ### 1. Build the MAMP database
 
 In a text file, save the following MAMP sequences (/MAMP_database/MAMP_elicitor_list.fasta):
 ```
 >csp22_consensus
AVGTVKWFNAEKGFGFITPDDG
>Elf18
SKEKFERTKPHVNVGTIG
>Flg22
QRLSTGSRINSAKDDAAGLQIA
>Nlp20
AIMYSWWFPKDSPVTGLGHR
```
This fasta file can be used to build a database to use blast to find if anything in the genome shares these sequences. To build the blast database, the below command was ran. Also this should be ran in the same folder as /MAMP_database/MAMP_elicitor_list.fasta. 

```
# to make blast db
makeblastdb -in MAMP_elicitor_list.fasta -parse_seqids -dbtype 'prot' -out MAMP_blast_db
````

We can then go through each protein fasta file and pull out the peptide from the annotation. In this case, we can use a bash loop to blast each file. 

```
for file in *.faa 
do echo "$file" 
blastp -task blastp-fast -query $file -db ./../../Mining_Known_MAMPs/MAMP_database/MAMP_blast_db -evalue 1e-4 -num_threads 4 -outfmt "6 qseqid sseqid pident evalue slen qstart qend length mismatch qseq" -out $file.txt 
done
```
```
for file in *.faa
do echo "$file"
blastp -task blastp-short -xdrop_gap_final 1000 -soft_masking false -query $file -db ./../../Mining_Known_MAMPs/MAMP_database/MAMP_blast_db -evalue 1e-4 -num_threads 4 -outfmt "6 qseqid sseqid pident evalue slen qstart qend length qseq" -out $file.txt
done
```

### 1. Build the MAMP database


 