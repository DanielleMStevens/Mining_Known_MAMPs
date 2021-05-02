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
  
  I have collected all the accession numbers as well as info about each one into two file stored in the Genome_accession_info directory. First we will need to install a package which will allow us to easily download the genomes by accession number from NCBI's refseq.
  
    conda install ncbi-genome-download
    
  Once that has installed, run the following command to download the genomes we are interested in. 

In order to get accessions for each key genus of bacteria

```bash
ncbi-genome-download -s refseq -g Agrobacterium --dry-run bacteria
```

```bash
ncbi-genome-download --assembly-accessions ./Genome_accession_info/Genome_accessions_to_download.txt -p 4 -r 2 -v --flat-output -F genbank,fasta,protein-fasta bacteria
```