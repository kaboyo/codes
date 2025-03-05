#http://ivory.idyll.org/blog/2023-snakemake-slithering-section-1.html
#snakemake for doing bioinformatics - a beginner's guide (part 1)
#Installation and setup!

```
#Setup and installation

#I suggest working in a new directory.
#You'll need to install snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) & 
# sourmash (https://sourmash.readthedocs.io/en/latest/#installing-sourmash). 
# We suggest using mamba, via miniforge/mambaforge (https://github.com/conda-forge/miniforge#mambaforge), for this.



```
conda create -c conda-forge -c bioconda -n snakemake snakemake
```
#Getting the data:

```
getting the data:
You'll need to download these three files:
GCF_000021665.1_ASM2166v1_genomic.fna.gz: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/021/665/GCF_000021665.1_ASM2166v1/GCF_000021665.1_ASM2166v1_genomic.fna.gz
GCF_000017325.1_ASM1732v1_genomic.fna.gz: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/325/GCF_000017325.1_ASM1732v1/GCF_000017325.1_ASM1732v1_genomic.fna.gz
GCF_000020225.1_ASM2022v1_genomic.fna.gz: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/225/GCF_000020225.1_ASM2022v1/GCF_000020225.1_ASM2022v1_genomic.fna.gz

```

#and rename them so that they are in a subdirectory genomes/ with the names:
```
GCF_000017325.1.fna.gz
GCF_000020225.1.fna.gz
GCF_000021665.1.fna.gz
```

#Creating a Snakefile for this project  

#Chapter 1 - snakemake runs programs for you!
#Here's a simple, useful snakemake workflow:
```
rule compare_genomes:
    message: "compare all input genomes using sourmash"
    shell: """
        sourmash sketch dna -p k=31 genomes/*.fna.gz --name-from-first 

        sourmash compare GCF_000021665.1.fna.gz.sig \
            GCF_000017325.1.fna.gz.sig GCF_000020225.1.fna.gz.sig \
            -o compare.mat

        sourmash plot compare.mat
    """
```

#Put it in a file called Snakefile, and run it with snakemake -j 1.
#This will produce the output file compare.mat.matrix.png which contains a similarity matrix and a dendrogram of the three genomes (see Figure 1).
*#copy the files to the local*
```
scp eo2r24@iridis5.soton.ac.uk:/scratch/eo2r24/SNAKEMAKE/tourmaline/00-data/metadata.tsv /Users/expeditoolimi/Documents/PRJNA358488/codes 
```