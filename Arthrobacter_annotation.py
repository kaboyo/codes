#good pictures: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.14036

#The workflows are available as open source software under the MIT license at https://github.com/pmenzel/ont-assembly-snake and https://github.com/pmenzel/score-assemblies.

Melange snakemake pipeline was adapted and used  (https://github.com/sandragodinhosilva/melange?tab=readme-ov-file)
#Melange is a genome annotation tool that enables the simultaneous annotation of large genome datasets using multiple databases: Pfam, COG, KEGG, CAZymes and/or MEROPS.
#Melange can handle unassembled and assembled sequencing data and amino acid sequences, with automatic download and configuration of necessary tools and databases.
#As a Snakemake pipeline, Melange is highly scalable, reproducible and has a transparent workflow, and can be used to annotate one to thousands of genomes, producing several easy-to-analyze tabular outputs.
#Melange derives its name from the French term "mélange", which signifies a blend or collection of diverse entities. This choice of name is a tribute to its ability to facilitate annotation using multiple databases.

For details, check: https://github.com/sandragodinhosilva/melange?tab=readme-ov-file


#A) Inputs
#Melange accepts 3 types of input files:

#A)unassembled (meta)genomic data (.fastq)
(meta)genome assemblies (.fna, .fasta, .ffn, .faa, .frn, .fa)
predicted amino acid sequences (.faa)
#If fastq files are inputted, Melange will convert them to fasta nucleotide files using the 
# EMBOSS tool seqret before annotation.

#B) Genome annotation
#B1) Gene calling and general annotation
When nucleotide files are submitted, Melange first performs a structural annotation step using Prokka v1.14.5 
[1] with default settings and outputs the corresponding translations into amino acid sequences in a fasta file. 
In addition to this output, which will be used in all subsequent steps, 
Prokka also generates other additional file formats, such as GenBank files, per genome.

#B2) Functional annotation
#Melange allows functional annotation of genomes with up to five databases: Pfam, COG, KEGG, CAZymes and MEROPS.

#Here are described the main characteristics of the annotation procedure with each database:

#Annotation databases

#Pfam: For the annotation with Pfam identifiers, a local database is created using HMMER v3.3 
# from the latest version of Pfam-A.hmm file (currently v35.0) downloaded from the downloaded 
# from the InterPro repository A local database is constructed using HMMER v3.3. Once the local database 
# has been created, query proteins are searched against it using the hmmscan function from the HMMER suite. 
# The best hit per ORF (cut-off: -E 1e-5) is selected. Selected references: [2,3].


#COG (Clusters of Orthologous Groups): The COG annotation procedure follows the cdd2cog v0.2 workflow. 
# First, several files are downloaded from the NCBI's FTP server, including a preformatted database of the 
# NCBI's Conserved Domain Database (CDD) COG distribution (2020 release). Query proteins are then blasted 
# against this database using reverse position-specific BLAST (rps-blast) function from the Blast+ v2.9.0 suite 
# and the results are parsed to a readable format with a Perl script (cdd2cog.pl). The best hit per ORF 
# (cut-off: -E 1e-5) is selected. Selected references: [4].

#KEGG (Kyoto Encyclopaedia of Genes and Genomes): To obtain the KEGG Orthology (KO) for protein identification, 
# the command line (CLI) version of KofamKoala(https://www.genome.jp/tools/kofamkoala/) - Kofamscan - is used. 
# Kofamscan performs K number assignments using hidden Markov model (HMM) profile search, which involves searching 
# query proteins against a customized HMM database of KOs (KEGG release 103.0). This database includes predefined 
# thresholds for individual KOs, resulting in more reliable assignments than sequence similarity searches. 
# Kofamscan uses the hmmsearch function from the HMMER suite to perform the search. Selected references: [5,6].


#CAZymes (Carbohydrate-active enzymes): The CAZymes annotation procedure uses the meta server dbCAN2, 
# specifically, the standalone version run_dbcan v2.0.11 implemented with default settings. 
# Run_dbcan is a tool that performs annotation of CAZymes using three different approaches: a HMMER v3.3 
# search against the dbCAN HMM database, a DIAMOND v0.9.32 search against the CAZy database, and the eCAMI algorithm. 
# For improved annotation accuracy, ORFs are only annotated with the respective CAZyme name if at least two database 
# searches were positive, as suggested by dbCAN2 authors in Zhang et al.. Selected references: [7,8].

#MEROPS: For the identification of ORFs encoding for peptidases and their inhibitors the "merops_scan.lib", 
# release 12.4 file is downloaded from MEROPS. Then makedblast is used to produce a local BLAST database. 
# Query aminoacid sequences are then searched for matches with this database with blastp. Selected references: [9].

#C) Outputs
#Melange produces several different output formats tailored to meet users' diverse needs, with almost no additional 
# computational cost. This is achieved by leveraging the output of each annotation database and transforming it 
# into different tables.


#In summary, three files with distinct data representation modes are created for each annotation type:

Counts
Presence/absence (PA)
Relative abundance
Outputs

#In these output tables, each row represents a database identifier (ID), and each column represents 
# an input (either nucleotide or amino acid (meta)genome files). While in counts (A), nij represents 
# the number of proteins or protein domains (depending on the database in use) identified with a certain
# ID for a given input, in the PA table (B), nij equal to 1 indicates the existence of a certain identifier 
# in the input, and 0 indicates its absence. In the relative abundance annotation table (C), nij represents 
# the normalized count of an ID per the total number of ORFs in each input.

#In addition to the annotation tables, Melange also provides:

#intermediate files - including different file types (e.g. GenBank, GFFF, etc);
#descriptive file (e.g. Pfam_description.csv) containing a summarized description of each annotation ID;
#statistics.csv - % of Orfs annotated with each database;
#folder Orf_per_genome: each genome has a unique file containing all orfs identified by Prokka and 
# the subsequent annotations with the selected functional databases;
#benchmark results - Melange automatically records running metrics using the Snakemake directive benchmark.


Usage
#This is a simple description on how to use melange. For more details, please see Melange documentation.

#Step 0: Install conda and Snakemake
#Conda and Snakemake are required to be able to use Melange.
#Conda is easy to install via its lightweight version Miniconda.
#After installing Conda, install Snakemake:

# As described in Snakemake documentation:
```
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

#Step 1: Clone Melange workflow

#To use Melange, you need a local copy of the workflow repository. Start by creating a clone of the repository:
```
git clone https://github.com/sandragodinhosilva/melange.git
```

#Step 2: Configure workflow
```
Configure the workflow according to your needs by editing the file config.yaml.
```
#Here you can select which databases (Pfam, COG, Kegg, CAZymes and/or MEROPS) are to be used.

#You can also define if input files are either fasta nucleotide files (e.g. fna, fa) or fasta aminoacid files.

#More information about configuration settings can be found at: config/README.md


#Step 3: Execute workflow
#Execute the workflow locally with N cores:
```
snakemake --use-conda --cores N
```

Execution on a cluster, example:
```
snakemake --use-conda --cluster qsub --jobs 8
```

#For more information about running on a computational cluster, please check snakemake documentation about it: https://snakemake.readthedocs.io/en/stable/executing/cluster.html

#Citing Melange
#At the moment, Melange does not have a publication describing its features (we are working on it). Please use a link to Melange Github when referring to this tool.

#Melange Contributions

#Sandra Godinho Silva 1,2 - MicroEcoEvo - iBB, IST.
#Tina Keller-Costa 1,2 - MicroEcoEvo - iBB, IST.
#Sanchita Kamath 3 - Microbial Systems Data Science - UFZ, Leipzig.
#Ulisses Nunes da Rocha 3 - Microbial Systems Data Science - UFZ, Leipzig.
#Rodrigo Costa 1,2 - MicroEcoEvo - iBB, IST.
#1 Institute for Bioengineering and Biosciences, Department of Bioengineering, Instituto Superior Técnico da Universidade de Lisboa, Lisbon, Portugal
#2 Associate Laboratory, Institute for Health and Bioeconomy, Instituto Superior Técnico, University of Lisbon, Lisbon, Portugal
#3 Department of Environmental Microbiology, Helmholtz Centre for Environmental Research – UFZ, Leipzig, Germany




#Project; Annotation of Arthrobacter genome

#Tools: Chatgpt search

################################################################
#Quality Check of FASTA Assembly
#Even though your sequences are already assembled, it’s important to check for errors before annotation.

```
 quast JDC6.genome.fa -o JDC6_quast_output
```
#OR
#CheckM – Evaluate genome completeness (if you have a reference lineage)
##Genome Size: Does it match expected size for your species?
#N50 & L50: Are contigs/scaffolds reasonable?
#GC Content: Is it consistent with reference genomes?
#Plasmid Circularization: Check if plasmid sequences are complete (use Bandage).
```
checkm lineage_wf -x fa JDC6.genome.fa.split/ checkm_output/ -t 8
```
#for checkm to work on the different genomes, we shall have to separate them main bacterial genome from the plasmid genome

#2. Split the Multi-Genome FASTA File (if needed)
#CheckM requires separate genome files. If all genomes are inside one single .fa file, split them using seqkit:

```
mkdir -p genomes
seqkit split -i -o genomes JDC6.genome.fa #adjust the names for the genomes
#This will create separate FASTA files in the genomes/ directory.
```

#If seqkit is not installed:
```
conda install -c bioconda seqkit
```

#Separate the Genome & Plasmid (if needed)
#If both the chromosome and plasmid(s) are in one FASTA file, split them manually or with a script.

#If they are separate already, skip this step.
#If they are in one FASTA file, identify the plasmid using:
```
plasmidfinder.py -i genome.fasta -o plasmidfinder_output
#https://gensoft.pasteur.fr/docs/platon/1.6/
#https://github.com/HubertTang/PLASMe
```


#Genome & Plasmid Annotation
#Use Prokka, a widely used annotation tool for bacterial genomes & plasmids.

#Annotate the genome
```
prokka --outdir genome_annotation --prefix genome JDC6.genome.fa.split/JDC6.genome.part_tig00000001.fa
```

#Annotate the plasmid
```
prokka --outdir plasmid_annotation --prefix plasmid JDC6.genome.fa.split/JDC6.genome.part_tig00000013.fa
```


#Functional Analysis (Optional, Post-Annotation)
#After annotation, you can analyze:

#Antimicrobial Resistance (AMR) Genes:

```
abricate --db resfinder genome.gff > amr_results.txt
#(Use CARD or ResFinder databases)
```

#Virulence Genes:

```
abricate --db vfdb genome.gff > virulence_results.txt
#Plasmid Replication & Mobilization Genes:
#Check plasmid-specific genes to confirm identity.
```
################################################################
conda create -n ncRNA -c bioconda barrnap tRNAscan-SE infernal crisprcasfinder aragorn
conda activate ncRNA



#Characterizing all non-coding RNAs (ncRNAs) in a bacterial genome involves several bioinformatics approaches and tools. Here’s a step-by-step guide:

#1. Assemble & Annotate the Genome
#Ensure that you have a high-quality genome assembly in FASTA format. If you haven't already assembled your genome, you can use tools like:

#SPAdes (for short-read assemblies)
#Canu (for long-read assemblies)
#Unicycler (for hybrid assemblies)
#2. Identify Non-Coding RNAs
#Several bioinformatics tools can be used for bacterial ncRNA prediction:

#a. rRNA (ribosomal RNA) Identification

#Barrnap – Rapid rRNA prediction
```
barrnap genome.fasta > rRNA.gff

#RNAmmer – More comprehensive rRNA annotation
```

#b. tRNA (transfer RNA) Prediction
tRNAscan-SE – The most widely used tool for tRNA identification
```
tRNAscan-SE -o tRNA_results.txt genome.fasta
```

#c. Small Regulatory RNAs (sRNAs) Prediction



#Infernal (Rfam Database) – Searches for known ncRNA families
```
cmsearch --cpu 4 --tblout sRNA_results.tbl Rfam.cm genome.fasta
sRNAPredict – Identifies novel bacterial sRNAs
#RNAz – Predicts thermodynamically stable ncRNAs
```

#d. CRISPR Arrays
#CRISPRCasFinder – Identifies CRISPR-Cas elements
```
CRISPRCasFinder -i genome.fasta -o crispr_results
```


##4. Visualization & Reporting
#Integrate results into a genome browser (e.g., IGV, Artemis)
#Use Python/R for custom visualization of ncRNA distribution
#There are many ways to visualize the genomes

#Here are the installation and usage commands for visualizing bacterial genomes using different tools.

#1. Circular Genome Visualization
(a) BRIG – Bacterial Genome Comparison
#Installation
```
# Download BRIG (Linux)
wget https://sourceforge.net/projects/brig/files/BRIG-0.95.zip
unzip BRIG-0.95.zip
cd BRIG-0.95

# Run BRIG
java -jar BRIG.jar
Usage
Load a reference genome (FASTA/GenBank)
Add other genomes to compare sequences
```
#Adjust settings for GC content, annotations, and identity thresholds
#(b) CGView – Circular Genome Visualization
#Installation & Usage (Java-based)
```
wget https://paulstothard.github.io/cgview/downloads/cgview.jar
java -jar cgview.jar
Load FASTA/GenBank files

java -jar cgview.jar -i genome.xml -o genome.png
#Customize tracks and labels
#Export the circular genome as SVG/PNG/PDF
```

#2. Linear Genome Visualization
#(a) Artemis – Genome Browser
#Installation
```
wget https://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.jar
java -jar artemis.jar
Usage
#Open the FASTA or GenBank genome file
##View CDS, GC content, and gene annotations
#Export annotated genome in EMBL or GenBank format
```

#(b) GenomeTools – Linear Genome Representation
#Installation
```
conda install -c bioconda genometools

# Convert GenBank file to graphical representation
gt sketch -format svg genome_output.svg genome.gbk
```

#(c) GenoPlotR – Genome Visualization in R
#Installation
#r
```
install.packages("devtools")
devtools::install_github("genouest/genoplotR")
Usage
r

library(genoPlotR)

# Load genome data
dna_seg1 <- read_dna_seg_from_file("genome1.gbk")
dna_seg2 <- read_dna_seg_from_file("genome2.gbk")

# Plot genomes
plot_gene_map(dna_segs=list(dna_seg1, dna_seg2))
```

#3. Comparative Genomics & Synteny Plots
#(a) Mauve – Multiple Genome Alignment
Installation
```
sudo apt install mauve
#Usage

mauve
#Load multiple bacterial genomes (FASTA format)
#Visualize genome rearrangements & conserved blocks
#Export alignment results as SVG/PNG

```

#(b) ACT – Artemis Comparison Tool
#Installation
```
wget https://ftp.sanger.ac.uk/pub/resources/software/artemis/act.jar
java -jar act.jar
Usage
Load two or more bacterial genomes
Compare gene synteny and structural variations
Export results for further analysis

```


#4. 3D & Web-Based Genome Visualization
#(a) JBrowse – Web-Based Genome Browser
#Installation
```
git clone https://github.com/GMOD/jbrowse.git
cd jbrowse
./setup.sh
Usage

./bin/flatfile-to-json.pl --gff genome_annotations.gff --trackLabel "Genes"
Open JBrowse in a browser and explore the genome interactively.
```
#(b) Circos – Advanced Circular Genome Visualization
#Installation
```
conda install -c bioconda circos
Usage

# Generate a circular genome plot
circos -conf circos.conf
#Customize Circos configuration (circos.conf) for visualizing GC content, SNPs, and synteny.
```
#star representation: https://github.com/robotoD/GenoVi/wiki/User-guide#keep-temporary-files-

#genovi -i genome.gbk -s draft -cs paradise --cogs_unclassified -bc yan51.cluster.local
#for the genome
#genovi -i genome.gbk -cs strong -s complete --size
#genovi -i genome.gbk -cs blossom -s draft --title 'filename'
#genovi -i input_test/Acinetobacter_radioresistens_DD78.gbff -cs par 
#genovi -i genome.gbk -cs paradise --scale linear --alignment
#genovi -i input_test/P_xenovorans_LB400.gbff -cs autumn -s draft -c
#genovi -i genome.gbk -cs autumn -s draft -c
#genovi -i genome.gbk -cs autumn -a A -s complete -c
#genovi -i genome.gbk -s complete --cogs met-
#genovi -i genome.gbk -s complete --cogs info-
#genovi -i genome.gbk -s complete --cogs cel-
#genovi -i genome.gbk -s complete --cogs poo-
#genovi -i genome.gbk -s complete --cogs poo-,met-,info-,cel-
#genovi -i genome.gbk -s complete --cogs QX
#genovi -i genome.gbk -s draft -t 'Arthrobacter sp.' --title_position top --italic_words 1 --size
#genovi -i genome.gbk -s complete -te -cs blossom

#plasmid
#genovi -i plasmid.gbk -s draft -cs paradise --cogs_unclassified -bc white
#genovi -i plasmid.gbk -cs strong -s complete --size
#genovi -i plasmid.gbk  -cs blossom -s draft --title 'filename'
#genovi -i plasmid.gbk  -cs paradise --scale linear --alignment '<' -s complete
#genovi -i plasmid.gbk  -cs autumn -s draft -c
#genovi -i plasmid.gbk -cs autumn -a A -s complete -c
#genovi -i plasmid.gbk -s complete --cogs met-
#genovi -i plasmid.gbk -s complete --cogs info-
#genovi -i plasmid.gbk -s complete --cogs cer-
#genovi -i plasmid.gbk -s complete --cogs poo-
#genovi -i plasmid.gbk -s complete --cogs QX
#genovi -i plasmid.gbk  -s complete --cogs 5
#genovi -i plasmid.gbk -s draft -t 'Arthrobacter sp.' --title_position top --italic_words 1 --size
#genovi -i plasmid.gbk  -s complete -te -cs blossom
#genovi -i plasmid.gbk -s complete -bc white-cs autumn -te --size
#genovi -i plasmid.gbk -s complete -bc white-cs autumn -te --size

################################################################

#Annotation Pipeline (Genome & Plasmid)
#Since the genome and plasmid are separate FASTA files, you should annotate them separately 
# to avoid misclassification of genes.

#
#✅ Bacterial genome FASTA (genome.fasta)
#✅ Plasmid FASTA (plasmid.fasta)

#Check if your FASTA headers are properly formatted. If needed, rename them:

```
sed -i 's/ .*//g' genome.fasta  # Remove spaces in headers
sed -i 's/ .*//g' plasmid.fasta
```

#2️⃣ Annotate Genome & Plasmid Separately
#Use Prokka (fast & widely used for bacterial annotation).

#Annotate Bacterial Genome:
```
prokka --outdir genome_annotation --prefix genome genome.fasta
```

#Annotate Plasmid:

```
prokka --outdir plasmid_annotation --prefix plasmid --plasmid plasmid.fasta
```
#✅ Output Files (.gff, .gbk, .faa, .ffn, etc.)

3#️⃣ Verify Annotations
#After annotation, check the functional elements:

```
grep -E "CDS|rRNA|tRNA" genome_annotation/genome.txt
```
#🔹 Check Plasmid-Specific Genes:

```
grep -E "replication|conjugation|mobility|tra" plasmid_annotation/plasmid.gff
```

#  🔹 Find Antimicrobial Resistance (AMR) Genes:

```
rgi main --input_sequence genome.fasta --output_file AMR_results.txt
```
#ResFinder (Plasmid AMR genes):
```
resfinder.py -i plasmid.fasta -o resfinder_output
```

#4️⃣ Comparative Genomics (Optional)
#If you want to compare your genome/plasmid with others:
    
```
#Pangenome Analysis: Roary, Panaroo
#Virulence Factor Analysis: VFDB
#Plasmid Typing: PlasmidFinder
#Final Output
#✔ genome.gff and plasmid.gff (Annotation results)
#✔ genome.gbk and plasmid.gbk (GenBank format for visualization)
#✔ Functional analysis (AMR genes, plasmid features, virulence factors)







################################################################
#Pangenome analysis
Used tools-PANGOLIN tools: https://github.com/labgem/PPanGGOLiN?tab=readme-ov-file
#check for other tools: https://github.com/UPHL-BioNGS/Grandeur/issues/188
# Panaroo: /https://gthlab.au/panaroo/#/gettingstarted/quickstart 
& 
#pirate: https://github.com/SionBayliss/PIRATE
################################################################

#Downloading the genomes from ncbi that are affilliated with Arthrobacter 

Details for all the priocess of download can be accessed here: https://github.com/kblin/ncbi-genome-download/

#Commad
```
ncbi-genome-download --genera "Arthrobacter" bacteria --parallel 2
```
#downloaded genomes are in gbff format
```
less genome.gbff
zcat genome.gbff.gz | less  # For compressed files
```
#the pangolin can implement the gbff runs, but i laboured to design a script to convert the genomes from gbff format to fasta format
```
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "biopython",
# ]
# ///
import os
import gzip
from Bio import SeqIO

# Define the parent directory to search in
SEARCH_DIR = "/scratch/eo2r24/Arthrobacter_pangenome/Genomes_arthrobacter/refseq"

# Define the destination directory for FASTA files
DEST_DIR = "/scratch/eo2r24/Arthrobacter_pangenome/Genomes_arthrobacter/Genomes-Arthrobacter-Fasta"

# Ensure the destination directory exists
os.makedirs(DEST_DIR, exist_ok=True)

# Function to convert .gbff.gz to .fasta
def convert_gbff_to_fasta(gbff_path, fasta_path):
    with gzip.open(gbff_path, "rt") as gbff_file:
        with open(fasta_path, "w") as fasta_file:
            for record in SeqIO.parse(gbff_file, "genbank"):
                SeqIO.write(record, fasta_file, "fasta")
    print(f"Converted: {gbff_path} _ {fasta_path}")

# Search for .gbff.gz files and convert them to FASTA
for root, _, files in os.walk(SEARCH_DIR):
    for file in files:
        if file.endswith(".gbff.gz"):
            gbff_path = os.path.join(root, file)

            # Define the FASTA file name and path
            fasta_filename = file.replace(".gbff.gz", ".fasta")
            fasta_path = os.path.join(DEST_DIR, fasta_filename)

            # Convert to FASTA
            convert_gbff_to_fasta(gbff_path, fasta_path)

#save the file as convert_gbff_to_fasta.py
#run: uv run convert_gbff_to_fasta.py
#if biopython is not installed, run: uv add --script convert_gbff_to_fasta.py biopython

```

#Run a complete pangenome analysis
Pangolin (https://github.com/labgem/PPanGGOLiN?tab=readme-ov-file) was used to construct the pagenomes 


#Installation

#PPanGGOLiN can be is easily installed via conda, accessible through the bioconda channel.
To ensure a smoother installation and avoid conflicting dependencies, it's highly recommended to create a dedicated environment for PPanGGOLiN:

# Install PPanGGOLiN into a new conda environment
```
conda create -n ppanggolin -c defaults -c conda-forge -c bioconda ppanggolin

# Check PPanGGOLiN install
conda activate ppanggolin
ppanggolin --version (version 2.2.1)
```


#list the files:
ls Genomes-Arthrobacter-Fasta/*.fasta >genome-files.txt
#Generate a list of files consdering sample and direction

import os

# Specify the directory where you want to search for FASTA files
directory = "/scratch/eo2r24/Arthrobacter_pangenome/Genomes_arthrobacter/Genomes-Arthrobacter-Fasta"  # Change this to your directory path

# Open a text file to write the output
with open("fasta_files_with_paths.txt", "w") as file:
    # Walk through the directory and subdirectories
    for dirpath, dirnames, filenames in os.walk(directory):
        # Filter for .fasta and .fa files
        for filename in filenames:
            if filename.endswith(('.fasta', '.fa')):
                # Extract the full path including extension
                full_path = os.path.join(dirpath, filename)
                # Remove the extension (.fasta or .fa) from the filename only
                filename_without_extension = os.path.splitext(filename)[0]
                # Write the filename (without extension) and the full path (with extension) to the text file
                file.write(f"{filename_without_extension}\t{full_path}\n")

print(f"List of FASTA files with paths (without extensions in filenames) saved to fasta_files_with_paths.txt")

#Design a slurm script to schedule the job

