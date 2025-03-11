
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

#details for all the priocess of download can be accessed here: https://github.com/kblin/ncbi-genome-download/
```
ncbi-genome-download --genera "Arthrobacter" bacteria --parallel 2
```


#. Basic Listing with Numbers
```
ls | nl
```

#Detailed Listing
```
ls -lh | nl
```

#Listing Including Hidden Files
```
ls -lah | nl
```
#Numbered Directory-Only Listing

```
ls -d */ | nl
```

Pangraphs: https://github.com/neherlab/pangraph
#Pangraph is currently undergoing a major migration between v0 and v1. 
# In this short transition period links and documentation may be inconsistent.

#Overview
##pangraph provides a command line interface to find homology amongst large collections of closely related genomes. 
#The core of the algorithm partitions each genome into blocks that represent a sequence interval related by vertical descent. 
# Each genome is then an ordered walk along blocks. The collection of all genomes form a graph that captures all 
# observed structural diversity. pangraph is a standalone tool useful to parsimoniously infer horizontal gene transfer events 
# within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many 
# related genomes.

#Keep the pipes on the hpc
module load tmux
tmux new -s mysession
Ctrl + B, then D #detach
tmux attach -t mysession reattacch
tmux attach-session -t mysession 

#https://github.com/swittouck/benchmark-scarap: Benchmark SCARAP

git clone https://github.com/SWittouck/benchmark-scarap.git
wget https://github.com/davidemms/Open_Orthobench/releases/download/v1.1/BENCHMARKS.tar.gz


```
compress genomes as gzip /*.fa

UNZIP commands: for zip in proteomes__*.zip ; do unzip $zip ; done
tar xzf *.tar.gz *.tar.gz
unzip xzf *.tar.gz *.tar.gz
zip xzf *.tar.gz *.tar.gz
```

#count comtigs:
#Using awk (Basic Contig Lengths)
```
awk '/^>/ {if (seqlen) print seqlen; print; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' JDC6.genome.fasta
```

#Using seqkit (Fast and Easy)
# If you have seqkit installed, use:

```
seqkit stats genome.fasta
```

#Using bioawk (Flexible Output)

```
bioawk -c fastx '{print $name, length($seq)}' genome.fasta
```

#Using awk (Basic Contig Lengths)

```
bioawk -c fastx '{print $name, length($seq)}' genome.fasta
```

#Using samtools faidx (Index and Summarize)
```
samtools faidx genome.fasta
#then display contig length
cut -f1,2 genome.fasta.fai
```

pipemake lr-reseq-assemble-hifiasm



################################################################

#pipeline for genome annotation

#1. Check Assembly Quality
#Run QUAST to assess genome quality:

```
quast -o quast_output JDC6_genome.fasta
```

#Number of contigs
#N50, L50 values
#GC content

#2. Identify Plasmids
#To check if your genome contains plasmids:

#Using PlasmidFinder
```
plasmidfinder.py -i JDC6_genome.fasta -o plasmid_output
```
#or Using BLAST (Against Plasmid Database)
```
blastn -query JDC6_genome.fasta -db plasmid_database -out plasmid_results.txt
```

#If plasmids are detected, separate them from the main genome.

#3. Annotate the Genome (Including Plasmids)

#Use Prokka for bacterial genome annotation:

```
prokka --outdir annotation_output --prefix JDC6_annotation JDC6_genome.fasta
```


#4. Verify Circularization (If Needed)
#Check if the plasmid is circular using Bandage:

```
bandage image assembly_graph.gfa output.png
```

snakemake -s Snakefile --cores 20 --use-conda
