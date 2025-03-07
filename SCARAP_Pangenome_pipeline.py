
#SCARAP: pangenome inference and comparative genomics of prokaryotes
#SCARAP is a toolkit with modules for various tasks related to comparative genomics of prokaryotes. 
# SCARAP has been designed to be fast and scalable. Its main feature is pangenome inference, 
# but it also has modules for direct core genome inference (without inferring the full pangenome), 
# subsampling representatives from a (large) set of genomes and constructing a concatenated core gene alignment ("supermatrix") 
# that can later be used for phylogeny inference. SCARAP has been designed for prokaryotes but
# should work for eukaryotic genomes as well. It can handle large genome datasets on a range of taxonomic levels;
# it has been tested on datasets with prokaryotic genomes from the species to the order level.

https://github.com/swittouck/scarap?tab=readme-ov-file



#conda installation
```
conda create -n scarap python=3.11
conda activate scarap
conda install bioconda::scarap
```
#pip installation
```
pip install scarap
```


#Inferring a pangenome
```
scarap pan ./faas ./pan -t 16
```

################################

#Obtaining data
#SCARAP works mainly with faa files: amino acid sequences of all (predicted) genes in a genome assembly. 
# You can obtain faa files in at least three ways:

#1You can run a gene prediction tool like Prodigal (https://github.com/hyattpd/Prodigal) on genome assemblies of your favorite strains, 
or a complete annotation pipeline such as Prokka (https://github.com/tseemann/prokka) or Bakta (https://github.com/oschwengers/bakta).

#2You can search your favorite taxon on NCBI genome (https://www.ncbi.nlm.nih.gov/datasets/genome/) and manually download assemblies in 
the following way: click on an assembly, click "Download", select "Protein (FASTA)" as file type and click "Download" again.

#3 Given a list of assembly accession numbers (i.e. starting with GCA/GCF), you can use ncbi-genome-download (https://github.com/kblin/ncbi-genome-download/) 
to download the corresponding faa files.
Given a list of accessions in a file called accessions.txt, you can use ncbi-genome-download to download faa files as follows:



for option 3, the ncbi-genome dowload was used to obtain the genimes affiliatted to the genus Arthrobacter
#NCBI Genome Downloading Scripts to get the genomes
#https://github.com/kblin/ncbi-genome-download/

#installation script
```
pip install ncbi-genome-download
```

#Alternatively, clone this repository from GitHub, then run (in a python virtual environment)

#installation script
```
pip install .
```

#If this fails on older versions of Python, try updating your pip tool first:

```
pip install --upgrade pip
```

and then rerun the ncbi-genome-download install.

#conda installation

#Alternatively, ncbi-genome-download is packaged in conda. 
# Refer the the Anaconda/miniconda site to install a distribution (highly recommended). With that installed one can do:
```
conda install -c bioconda ncbi-genome-download
```  