
**#use the snakemake script defined basing on https://github.com/aomlomics/tourmaline?tab=readme-ov-file (Tourmaline pipeline)**
*this pipeline provides seamless analysis of amplicon read data, using the QIIME2 and the associated wrapped packages— the pipeline manages commands, input and output using the snakemake workflow management system. the current implementation uses/suport qiime2-2023.5 (https://docs.qiime2.org/2023.5/).*

#connect to the workstation 
```
ssh xxxxx@xxxxxxx
```

**#interactive runs on work station:**
```
sinteractive --partition=highmem --cpus-per-task=20 --mem=8000 --time=01:03:04
```


**#To run Tourmaline natively on a Mac (Intel) or Linux system, start with a Conda installation of Snakemake (https://snakemake.readthedocs.io/en/stable/)**

```
conda create -c conda-forge -c bioconda -n snakemake snakemake-minimal
```


**#Install qiime2-2023.5**
```
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-osx-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-osx-conda.yml
```
**#Activate the qiime2-2023.5 environment and install the other Conda- or PIP-installable dependencies:**

```
conda activate qiime2-2023.5
conda install -c conda-forge -c bioconda biopython muscle clustalo tabulate
conda install -c conda-forge deicode
pip install empress
qiime dev refresh-cache
conda install -c bioconda bioconductor-msa bioconductor-odseq
```


**#Apple Silicon Macs**
**#Follow these instructions for Macs with M1/M2 chips**

**#First, set your Terminal application to run in Rosetta mode.**

```
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-osx-conda.yml
CONDA_SUBDIR=osx-64 conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-osx-conda.yml
conda activate qiime2-2023.5
conda config --env --set subdir osx-64
```
***#Docker installation, please see the documentation:https://github.com/aomlomics/tourmaline?tab=readme-ov-file********************************


***SETUP:********************************
*#If this is your first time running Tourmaline, you'll need to set up your directory. Simplified instructions are below, but see the Wiki's Setup (https://github.com/aomlomics/tourmaline/wiki/3-Setup) page for complete instructions.*


```
#How to create an env from the previous one (cloning)
git clone https://github.com/aomlomics/tourmaline.git
```

*#Setup for the test data*
*#The test data (16 samples of paired-end 16S rRNA data with 1000 sequences per sample) comes with your cloned copy of #Tourmaline. It's fast to run and will verify that you can run the workflow.*

**#Download reference database sequence and taxonomy files, named refseqs.qza and reftax.qza (QIIME 2 archives), in 01-imported:********************************
*#make sure the folder named imported is in the tourmarine directory*


```
cd tourmaline/01-imported
wget https://data.qiime2.org/2023.5/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2023.5/common/silva-138-99-tax-515-806.qza
ln -s silva-138-99-seqs-515-806.qza refseqs.qza
ln -s silva-138-99-tax-515-806.qza reftax.qza
```

**#Edit FASTQ manifests manifest_se.csv and manifest_pe.csv in 00-data so file paths match the location of your tourmalinedirectory. In the command below, replace /path/to with the location of your tourmaline directory—or skip this step if you are using the Docker container and you cloned tourmaline into /data:********************************
#as a by the way, just in case the file paths are not well-formed
```
cd ../00-data
cat manifest_pe.csv | sed 's|/data/tourmaline|/path/to/tourmaline|' > temp && mv temp manifest_pe.csv 
cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
```


**#Prepare the sequence reads (manifest file)********************************
################################################################
#code to organise the reads into a manifest file (SE reads)
```
#!/usr/bin/env python
import os
import glob
import pandas as pd
import argparse

def generate_manifest(input_dir, output_dir):
    # Ensure the output directory exists, create if it doesn't
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Search for forward read files (_R1.fq.gz)
    forward_reads = sorted(glob.glob(os.path.join(input_dir, "*_R1.fq.gz")))

    # Check if any forward reads exist
    if not forward_reads:
        print("Error: No forward reads (_R1) found in the specified directory.")
        exit(1)

    # Prepare data for manifest file
    manifest_data = []

    # Add forward reads to the manifest data
    for file in forward_reads:
        sample_name = os.path.basename(file).replace("_R1.fq.gz", "")
        absolute_path = os.path.abspath(file)
        manifest_data.append([sample_name, absolute_path, "forward"])

    # Convert to DataFrame and save as CSV with commas
    manifest_df = pd.DataFrame(manifest_data, columns=["sample-id", "absolute-filepath", "direction"])

    # Define the output file path
    output_file = os.path.join(output_dir, "manifest.tsv")

    # Write the manifest to the output directory
    manifest_df.to_csv(output_file, sep=",", index=False)

    print(f"Manifest file 'manifest.tsv' has been successfully generated in {output_dir}.")

if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a manifest file for single-end reads (forward reads).")
    
    # Input directory argument
    parser.add_argument('input_dir', type=str, help="Directory containing the FASTQ files.")
    
    # Output directory argument
    parser.add_argument('output_dir', type=str, help="Directory to save the output manifest file.")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to generate the manifest
    generate_manifest(args.input_dir, args.output_dir)

```
*# commands*
python prepare_SE_table.py /path/to/raw_reads /path/to/output_folder
 or
 ./prepare_SE_table.py /path/to/raw_reads /path/to/output_folder




**#consider the PE reads********************************
```
#!/usr/bin/env python
import os
import glob
import pandas as pd
import argparse

def generate_manifest(input_dir, output_dir):
    # Ensure the output directory exists, create if it doesn't
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Search for forward (_R1.fastq.gz) and reverse (_R2.fastq.gz) reads using wildcards
    forward_reads = sorted(glob.glob(os.path.join(input_dir, "*_R1.fastq.gz")))
    reverse_reads = sorted(glob.glob(os.path.join(input_dir, "*_R2.fastq.gz")))

    # Check if any reads exist
    if not forward_reads and not reverse_reads:
        print("Error: No FASTQ files found in the specified directory.")
        exit(1)

    # Prepare data for manifest file
    manifest_data = []

    # Add forward reads
    for file in forward_reads:
        sample_name = os.path.basename(file).replace("_R1.fastq.gz", "")
        absolute_path = os.path.abspath(file)
        manifest_data.append([sample_name, absolute_path, "forward"])

    # Add reverse reads
    for file in reverse_reads:
        sample_name = os.path.basename(file).replace("_R2.fastq.gz", "")
        absolute_path = os.path.abspath(file)
        manifest_data.append([sample_name, absolute_path, "reverse"])

    # Convert to DataFrame and save as CSV with commas
    manifest_df = pd.DataFrame(manifest_data, columns=["sample-id", "absolute-filepath", "direction"])

    # Define the output file path
    output_file = os.path.join(output_dir, "manifest.tsv")

    # Write the manifest to the output directory
    manifest_df.to_csv(output_file, sep=",", index=False)

    print(f"Manifest file 'manifest.tsv' has been successfully generated in {output_dir}.")


if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a manifest file for amplicon sequencing analysis.")
    
    # Input directory argument
    parser.add_argument('input_dir', type=str, help="Directory containing the FASTQ files.")
    
    # Output directory argument
    parser.add_argument('output_dir', type=str, help="Directory to save the output manifest file.")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to generate the manifest
    generate_manifest(args.input_dir, args.output_dir)

```
**#command**********************************
```
python generate_manifest.py /path/to/raw_reads /path/to/output_folder
or
python generate_manifest.py /path/to/raw_reads /path/to/output_folder
```


**#prepare the samples as described in the here** 
**#Now edit, replace, or store the required input files as described here:**


```
#1.Edit or replace the metadata file 00-data/metadata.tsv. The first column header should be "sample_name", with sample names matching the FASTQ manifest(s), and additional columns containing any relevant metadata for your samples. You can use a spreadsheet editor like Microsoft Excel or LibreOffice, but make sure to export the output in tab-delimited text format (.tsv file).

#2.Prepare FASTQ data:

    #Option 1: Edit or replace the FASTQ manifests 00-data/manifest_pe.csv (paired-end) and/or 00-data/manifest_se.csv (single-end). Ensure that (1) file paths in the column "absolute-filepath" point to your .fastq.gz files (they can be anywhere on your computer) and (2) sample names match the metadata file. You can use a text editor such as Sublime Text, nano, vim, etc.

    *#Option 2: Store your pre-imported FASTQ .qza files as 01-imported/fastq_pe.qza (paired-end) and/or 01-imported/fastq_se.qza (single-end).

#3.Prepare reference database:
    #Option 1: Store the reference FASTA and taxonomy files as 00-data/refseqs.fna and 00-data/reftax.tsv.

    #Option 2: Store the pre-imported reference FASTA and taxonomy .qza files as 01-imported/refseqs.qza and 01-imported/reftax.qza.

#4. Edit the configuration file config.yaml to set DADA2 and/or Deblur parameters (sequence truncation/trimming, sample pooling, chimera removal, etc.), rarefaction depth, taxonomic classification method, and other parameters. This YAML (yet another markup language) file is a regular text file that can be edited in Sublime Text, nano, vim, etc.
#5. Go to Run Snakemake.
```
################################

#Run Snakemake

**#Tourmaline is now run within the snakemake conda environment, not the qiime2-2023.5 environment.**

```
conda activate snakemake
```


*#Shown here is the DADA2 paired-end workflow. See the Wiki's Run (https://github.com/aomlomics/tourmaline/wiki/4-Run) page for complete instructions on all steps, denoising methods, and filtering modes.*

Note that any of the commands below can be run with various options, 
including 
```
--printshellcmds 
```
*#to see the shell commands being executed* and 
```
--dryrun 
```
*#to display which rules would* be run but not execute them. 

**#To generate a graph of the rules that will be run from any Snakemake command, see the section "Directed acyclic graph (DAG)" on the Run (https://github.com/aomlomics/tourmaline/wiki/4-Run) page.**

#Always include the 

```
--use-conda
```
*#option*

*#From the tourmaline directory (which you may rename), run Snakemake with the denoise rule as the target, changing the number of cores to match your machine:*

```
#test
snakemake --use-conda dada2_pe_denoise --cores 4 --dryrun 

#test-datab
snakemake --use-conda dada2_pe_denoise --cores 4 --printshellcmds

```
*#Pausing after the denoise step allows you to make changes before proceeding:*

*#Check the table summaries and representative sequence lengths to determine if DADA2 or Deblur parameters need to be modified. If so, you can rename or delete the output directories and then rerun the denoise rule.*

*#View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters "alpha_max_depth" and "core_sampling_depth" in config.yaml.*

*#Decide whether to filter your feature table and representative sequences by taxonomy or feature ID. After the taxonomy step, you can examine the taxonomy summary and bar plot to aid your decision. If you do filter your data, all output from that point on will go in a separate folder so you can compare output with and without filtering.*

**#Unfiltered mode**

*#Continue the workflow without filtering (for now). If you are satisfied with your parameters and files, run the taxonomy rule (for unfiltered data):*

```
snakemake --use-conda dada2_pe_taxonomy_unfiltered --cores 4 --printshellcmds

```
*#Next, run the diversity rule (for unfiltered data):*

```
snakemake --use-conda dada2_pe_diversity_unfiltered --cores 4 --printshellcmds
```

#Finally, run the report rule (for unfiltered data):

```
snakemake --use-conda dada2_pe_report_unfiltered --cores 4 --printshellcmds
```


*#Filtered mode*

*#After viewing the unfiltered results—the taxonomy summary and taxa barplot, the representative sequence summary plot and table, or the list of unassigned and potential outlier representative sequences—the user may wish to filter (remove) certain taxonomic groups or representative sequences.***

*#If so, the user should first check the following parameters and/or files:*

    #copy2-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv to 00-data/repseqs_to_filter_dada2-pe.tsv to filter outliers, or manually include feature IDs in 00-data/repseqs_to_filter_dada2-pe.tsv to filter those feature IDs (change "dada2-pe" to "dada2-se" or "deblur-se" as appropriate);

    #exclude_terms in config.yaml – add taxa to exclude from representative sequences, if desired;

    #repseq_min_length and repseq_max_length in config.yaml – set minimum and/or maximum lengths for filtering representative sequences, if desired;

    #repseq_min_abundance and repseq_min_prevalence in config.yaml – set minimum abundance and/or prevalence values for filtering representative sequences, if desired.

*#Now we are ready to filter the representative sequences and feature table, generate new summaries, and generate a new taxonomy bar plot, by running the taxonomy rule (for filtered data):*


```
snakemake --use-conda dada2_pe_taxonomy_filtered --cores 4  --printshellcmds
```

**##From the tourmaline directory (which you may rename), run Snakemake with the denoise rule as the target, changing the number of cores to match your machine:*******************************
