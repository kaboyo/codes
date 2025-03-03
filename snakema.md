
#The tourmaline pipeline for ampliseq analysis using snakemake snakemake 

#site:https://github.com/aomlomics/tourmaline


ssh commad ssh xxxxx@iridis5.soton.ac.uk 
password:


copy files between directoris: scp -r /Users/expeditoolimi/Documents/NEXTFLOW/PRJNA358488 eo2r24@iridis5.soton.ac.uk:/scratch/eo2r24/nexflow/PRJNA358488 



#snakemake
Follow the procedures here....

https://github.com/aomlomics/tourmaline

&https://github.com/olabiyi/snakemake-workflow-qiime2?tab=readme-ov-file (good)

#Tourmaline
#Tourmaline is an amplicon sequence processing workflow for Illumina sequence data that uses QIIME 2 and the software packages it wraps. Tourmaline manages commands, inputs, and outputs using the Snakemake workflow management system.

The current version of Tourmaline supports qiime2-2023.5. To use previous versions of Qiime2, check out previous Tourmaline versions under Releases.

If you would like to test the development version of Tourmaline 2, check out the develop branch of this repository!

#natively install snakemake
```
conda create -c conda-forge -c bioconda -n snakemake snakemake-minimal
```

#Then install QIIME 2 with conda (for Linux, change "osx" to "linux"):
```
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-osx-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-osx-conda.yml
```
#Activate the qiime2-2023.5 environment and install the other Conda- or PIP-installable dependencies:
```
conda activate qiime2-2023.5
conda install -c conda-forge -c bioconda biopython muscle clustalo tabulate
conda install -c conda-forge deicode
pip install empress
qiime dev refresh-cache
conda install -c bioconda bioconductor-msa bioconductor-odseq
```
#set up
#If this is your first time running Tourmaline, you'll need to set up your directory. Simplified instructions are below, but see the Wiki's Setup page for complete instructions.

```
git clone https://github.com/aomlomics/tourmaline.git
```

#set up test data sets


```
#make a directory 01-imported inside the tourmaline folder
cd tourmaline/01-imported
wget https://data.qiime2.org/2023.5/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2023.5/common/silva-138-99-tax-515-806.qza
ln -s silva-138-99-seqs-515-806.qza refseqs.qza
ln -s silva-138-99-tax-515-806.qza reftax.qza
```

#Edit FASTQ manifests manifest_se.csv and manifest_pe.csv in 00-data so file paths match the location of your tourmaline directory. In the command below, replace /path/to with the location of your tourmaline directory—or skip this step if you are using the Docker container and you cloned tourmaline into /data:

```
cd ../00-data
cat manifest_pe.csv | sed 's|/tourmaline|/path/to/tourmaline|' > temp && mv temp manifest_pe.csv 
cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
```

Meaning of the code: cat manifest_pe.csv | sed 's|/tourmaline|/path/to/tourmaline|' > temp && mv temp manifest_pe.csv 

```
#This command sequence is used to modify the contents of the file manifest_pe.csv. Here's a breakdown of what each part does:
#cat manifest_pe.csv: This command outputs the contents of the file manifest_pe.csv.

#|: This pipe symbol takes the output from the cat command and passes it as input to the next command.
#sed 's|/data/tourmaline|/path/to/tourmaline|': This sed command performs a substitution. It replaces all occurrences of the string /data/tourmaline with /path/to/tourmaline in the input it receives.

#> temp: This redirects the output of the sed command to a temporary file named temp.

#&&: This logical AND operator ensures that the next command (mv temp manifest_pe.csv) is only executed if the previous command (sed 's|/data/tourmaline|/path/to/tourmaline|' > temp) succeeds.

#mv temp manifest_pe.csv: This command renames (moves) the temporary file temp to manifest_pe.csv, effectively overwriting the original manifest_pe.csv with the modified content.

Meaning of the code: cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
#This command sequence is used to filter the contents of the file manifest_pe.csv and create a new file manifest_se.csv without lines containing the word "reverse". Here's a breakdown of what each part does:
#cat manifest_pe.csv: This command outputs the contents of the file manifest_pe.csv.
#|: This pipe symbol takes the output from the cat command and passes it as input to the next command.

#grep -v "reverse": This grep command searches for lines that do not contain the word "reverse". The -v option inverts the match, so it selects lines that do not match the specified pattern.

#> manifest_se.csv: This redirects the output of the grep command to a new file named manifest_se.csv.
```


Go to Run Snakemake.
```
#Setup for your data
#Before setting up to run your own data, please note:

#Symbolic links can be used for any of the input files, which may be useful for large files (e.g., the FASTQ and reference database .qza files).
#If you plan on using Deblur, sample names must not contain underscores (only alphanumerics, dashes, and/or periods).
#Now edit, replace, or store the required input files as described here:

#Edit or replace the metadata file 00-data/metadata.tsv. The first column header should be "sample_name", with sample names matching the FASTQ manifest(s), and additional columns containing any relevant metadata for your samples. You can use a spreadsheet editor like Microsoft Excel or LibreOffice, but make sure to export the output in tab-delimited text format.

#Prepare FASTQ data:
#Option 1: Edit or replace the FASTQ manifests 00-data/manifest_pe.csv (paired-end) and/or 00-data/manifest_se.csv (single-end). Ensure that (1) file paths in the column "absolute-filepath" point to your .fastq.gz files (they can be anywhere on your computer) and (2) sample names match the metadata file. You can use a text editor such as Sublime Text, nano, vim, etc.
#Option 2: Store your pre-imported FASTQ .qza files as 01-imported/fastq_pe.qza (paired-end) and/or 01-imported/fastq_se.qza (single-end).

#Prepare reference database:
#Option 1: Store the reference FASTA and taxonomy files as 00-data/refseqs.fna and 00-data/reftax.tsv.
#Option 2: Store the pre-imported reference FASTA and taxonomy .qza files as 01-imported/refseqs.qza and 01-imported/reftax.qza.

#Edit the configuration file config.yaml to set DADA2 and/or Deblur parameters (sequence truncation/trimming, sample pooling, chimera removal, etc.), rarefaction depth, taxonomic classification method, and other parameters. This YAML (yet another markup language) file is a regular text file that can be edited in Sublime Text, nano, vim, etc.
Go to Run Snakemake.
Run Snakemake

```



Tourmaline is now run within the snakemake conda environment, not the qiime2-2023.5 environment.
```
conda activate snakemake

```


#Shown here is the DADA2 paired-end workflow. See the Wiki's Run (https://github.com/aomlomics/tourmaline/wiki/4-Run) page for complete instructions on all steps, denoising methods, and filtering modes.

#Note that any of the commands below can be run with various options, including --printshellcmds to see the shell commands being executed and --dryrun to display which rules would be run but not execute them. To generate a graph of the rules that will be run from any Snakemake command, see the section "Directed acyclic graph (DAG)" on the Run page. Always include the --use-conda option.

From the tourmaline directory (which you may rename), run Snakemake with the denoise rule as the target, changing the number of cores to match your machine:
```
#snakemake --use-conda dada2_pe_denoise --cores 4
snakemake --use-conda dada2_pse_denoise --cores 4

#change cores when you like
```


#Pausing after the denoise step allows you to make changes before proceeding:

#Check the table summaries and representative sequence lengths to determine if DADA2 or Deblur parameters need to be modified. If so, you can rename or delete the output directories and then rerun the denoise rule.

#View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters "alpha_max_depth" and "core_sampling_depth" in config.yaml.

#Decide whether to filter your feature table and representative sequences by taxonomy or feature ID. After the taxonomy step, you can examine the taxonomy summary and bar plot to aid your decision. If you do filter your data, all output from that point on will go in a separate folder 

#Wso you can compare output with and without filtering.
#Unfiltered mode
#Continue the workflow without filtering (for now). If you are satisfied with your parameters and files, run the taxonomy rule (for unfiltered data):

```
snakemake --use-conda dada2_pe_taxonomy_unfiltered --cores 4
```

#Next, run the diversity rule (for unfiltered data):

```
snakemake --use-conda dada2_pe_diversity_unfiltered --cores 4
```

Finally, run the report rule (for unfiltered data):

```
snakemake --use-conda dada2_pe_report_unfiltered --cores 4

```





#Filtered mode
#After viewing the unfiltered results—the taxonomy summary and taxa barplot, the representative sequence summary plot and table, or the list of unassigned and potential outlier representative sequences—the user may wish to filter (remove) certain taxonomic groups or representative sequences. If so, the user should first check the following parameters and/or files:

#copy 2-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv to 00-data/repseqs_to_filter_dada2-pe.tsv to filter outliers, or manually include feature IDs in 00-data/repseqs_to_filter_dada2-pe.tsv to filter those feature IDs (change "dada2-pe" to "dada2-se" or "deblur-se" as appropriate);

#exclude_terms in config.yaml – add taxa to exclude from representative sequences, if desired;
repseq_min_length and repseq_max_length in config.yaml – set minimum and/or maximum lengths for filtering representative sequences, if desired;

#repseq_min_abundance and repseq_min_prevalence in config.yaml – set minimum abundance and/or prevalence values for filtering representative sequences, if desired.

#Now we are ready to filter the representative sequences and feature table, generate new summaries, and generate a new taxonomy bar plot, by running the taxonomy rule (for filtered data):

```
snakemake --use-conda dada2_pe_taxonomy_filtered --cores 4
```

Next, run the diversity rule (for filtered data):


```
snakemake --use-conda dada2_pe_diversity_filtered --cores 4

```

#Finally, run the report rule (for filtered data):

```
snakemake --use-conda dada2_pe_report_filtered --cores 1

```


```
View output
View report and output files
Open your HTML report (e.g., 03-reports/report_dada2-pe_unfiltered.html) in Chrome{target="_blank"} or Firefox{target="_blank"}. To view the linked files:

QZV (QIIME 2 visualization): click to download, then drag and drop in https://view.qiime2.org{target="_blank"}. Empress trees (e.g., rooted_tree.qzv) may take more than 10 minutes to load.
TSV (tab-separated values): click to download, then open in Microsoft Excel or Tabview (command line tool that comes with Tourmaline).
PDF (portable document format): click to open and view in new tab.
Downloaded files can be deleted after viewing because they are already stored in your Tourmaline directory.

More tips
Troubleshooting

#The whole workflow with test data should take ~3–5 minutes to complete. A normal dataset may take several hours to complete.

#If any of the above commands don't work, read the error messages carefully, try to figure out what went wrong, and attempt to fix the offending file. A common issue is the file paths in your FASTQ manifest file need to be updated.


#If you are running in a Docker container and you get an error like "Signals.SIGKILL: 9", you probably need to give Docker more memory. See the Wiki section on Installation options (https://github.com/aomlomics/tourmaline#:~:text=View%20output,the%20command%20below%3A).

Power tips

#The whole workflow can be run with just the command snakemake dada2_pe_report_unfiltered (without filtering representative sequences) or snakemake dada2_pe_report_filtered (after filtering representative sequences). Warning: If your parameters are not optimized, the results will be suboptimal (garbage in, garbage out).

#If you want to make a fresh run and not save the previous output, simply delete the output directories (e.g., 02-output-{method}-{filter} and 03-report) generated in the previous run. If you want to save these outputs and rerun with different parameters, you can change the name of the output directories and report files to something informative and leave them in the Tourmaline directory.

#You can always delete any file you want to regenerate. Then there are several ways to regenerate it: run snakemake FILE and Snakemake will determine which rules (commands) need to be run to generate that file; or, run snakemake RULE where the rule generates the desired file as output.

#If you've run Tourmaline on your dataset before, you can speed up the setup process and initialize a new Tourmaline directory (e.g., tourmaline-new) with the some of the files and symlinks of the existing one (e.g., tourmaline-existing) using the command below:
```
```
cd /path/to/tourmaline-new
scripts/initialize_dir_from_existing_tourmaline_dir.sh /path/to/tourmaline-existing
```



#The tourmaline pipeline for ampliseq analysis using snakemake snakemake 

#site:https://github.com/aomlomics/tourmaline


ssh commad ssh xxxxx@iridis5.soton.ac.uk 
password:


copy files between directoris: scp -r /Users/expeditoolimi/Documents/NEXTFLOW/PRJNA358488 eo2r24@iridis5.soton.ac.uk:/scratch/eo2r24/nexflow/PRJNA358488 



#snakemake
Follow the procedures here....

https://github.com/aomlomics/tourmaline

&https://github.com/olabiyi/snakemake-workflow-qiime2?tab=readme-ov-file (good)

#Tourmaline
#Tourmaline is an amplicon sequence processing workflow for Illumina sequence data that uses QIIME 2 and the software packages it wraps. Tourmaline manages commands, inputs, and outputs using the Snakemake workflow management system.

The current version of Tourmaline supports qiime2-2023.5. To use previous versions of Qiime2, check out previous Tourmaline versions under Releases.

If you would like to test the development version of Tourmaline 2, check out the develop branch of this repository!

#natively install snakemake
```
conda create -c conda-forge -c bioconda -n snakemake snakemake-minimal
```

#Then install QIIME 2 with conda (for Linux, change "osx" to "linux"):
```
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-osx-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-osx-conda.yml
```
#Activate the qiime2-2023.5 environment and install the other Conda- or PIP-installable dependencies:
```
conda activate qiime2-2023.5
conda install -c conda-forge -c bioconda biopython muscle clustalo tabulate
conda install -c conda-forge deicode
pip install empress
qiime dev refresh-cache
conda install -c bioconda bioconductor-msa bioconductor-odseq
```
#set up
#If this is your first time running Tourmaline, you'll need to set up your directory. Simplified instructions are below, but see the Wiki's Setup page for complete instructions.

```
git clone https://github.com/aomlomics/tourmaline.git
```

#set up test data sets


```
#make a directory 01-imported inside the tourmaline folder
cd tourmaline/01-imported
wget https://data.qiime2.org/2023.5/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2023.5/common/silva-138-99-tax-515-806.qza
ln -s silva-138-99-seqs-515-806.qza refseqs.qza
ln -s silva-138-99-tax-515-806.qza reftax.qza
```

#Edit FASTQ manifests manifest_se.csv and manifest_pe.csv in 00-data so file paths match the location of your tourmaline directory. In the command below, replace /path/to with the location of your tourmaline directory—or skip this step if you are using the Docker container and you cloned tourmaline into /data:

```
cd ../00-data
cat manifest_pe.csv | sed 's|/tourmaline|/path/to/tourmaline|' > temp && mv temp manifest_pe.csv 
cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
```

Meaning of the code: cat manifest_pe.csv | sed 's|/tourmaline|/path/to/tourmaline|' > temp && mv temp manifest_pe.csv 

```
#This command sequence is used to modify the contents of the file manifest_pe.csv. Here's a breakdown of what each part does:
#cat manifest_pe.csv: This command outputs the contents of the file manifest_pe.csv.

#|: This pipe symbol takes the output from the cat command and passes it as input to the next command.
#sed 's|/data/tourmaline|/path/to/tourmaline|': This sed command performs a substitution. It replaces all occurrences of the string /data/tourmaline with /path/to/tourmaline in the input it receives.

#> temp: This redirects the output of the sed command to a temporary file named temp.

#&&: This logical AND operator ensures that the next command (mv temp manifest_pe.csv) is only executed if the previous command (sed 's|/data/tourmaline|/path/to/tourmaline|' > temp) succeeds.

#mv temp manifest_pe.csv: This command renames (moves) the temporary file temp to manifest_pe.csv, effectively overwriting the original manifest_pe.csv with the modified content.

Meaning of the code: cat manifest_pe.csv | grep -v "reverse" > manifest_se.csv
#This command sequence is used to filter the contents of the file manifest_pe.csv and create a new file manifest_se.csv without lines containing the word "reverse". Here's a breakdown of what each part does:
#cat manifest_pe.csv: This command outputs the contents of the file manifest_pe.csv.
#|: This pipe symbol takes the output from the cat command and passes it as input to the next command.

#grep -v "reverse": This grep command searches for lines that do not contain the word "reverse". The -v option inverts the match, so it selects lines that do not match the specified pattern.

#> manifest_se.csv: This redirects the output of the grep command to a new file named manifest_se.csv.
```


Go to Run Snakemake.
```
#Setup for your data
#Before setting up to run your own data, please note:

#Symbolic links can be used for any of the input files, which may be useful for large files (e.g., the FASTQ and reference database .qza files).
#If you plan on using Deblur, sample names must not contain underscores (only alphanumerics, dashes, and/or periods).
#Now edit, replace, or store the required input files as described here:

#Edit or replace the metadata file 00-data/metadata.tsv. The first column header should be "sample_name", with sample names matching the FASTQ manifest(s), and additional columns containing any relevant metadata for your samples. You can use a spreadsheet editor like Microsoft Excel or LibreOffice, but make sure to export the output in tab-delimited text format.

#Prepare FASTQ data:
#Option 1: Edit or replace the FASTQ manifests 00-data/manifest_pe.csv (paired-end) and/or 00-data/manifest_se.csv (single-end). Ensure that (1) file paths in the column "absolute-filepath" point to your .fastq.gz files (they can be anywhere on your computer) and (2) sample names match the metadata file. You can use a text editor such as Sublime Text, nano, vim, etc.
#Option 2: Store your pre-imported FASTQ .qza files as 01-imported/fastq_pe.qza (paired-end) and/or 01-imported/fastq_se.qza (single-end).

#Prepare reference database:
#Option 1: Store the reference FASTA and taxonomy files as 00-data/refseqs.fna and 00-data/reftax.tsv.
#Option 2: Store the pre-imported reference FASTA and taxonomy .qza files as 01-imported/refseqs.qza and 01-imported/reftax.qza.

#Edit the configuration file config.yaml to set DADA2 and/or Deblur parameters (sequence truncation/trimming, sample pooling, chimera removal, etc.), rarefaction depth, taxonomic classification method, and other parameters. This YAML (yet another markup language) file is a regular text file that can be edited in Sublime Text, nano, vim, etc.
Go to Run Snakemake.
Run Snakemake

```



Tourmaline is now run within the snakemake conda environment, not the qiime2-2023.5 environment.
```
conda activate snakemake

```


#Shown here is the DADA2 paired-end workflow. See the Wiki's Run (https://github.com/aomlomics/tourmaline/wiki/4-Run) page for complete instructions on all steps, denoising methods, and filtering modes.

#Note that any of the commands below can be run with various options, including --printshellcmds to see the shell commands being executed and --dryrun to display which rules would be run but not execute them. To generate a graph of the rules that will be run from any Snakemake command, see the section "Directed acyclic graph (DAG)" on the Run page. Always include the --use-conda option.

From the tourmaline directory (which you may rename), run Snakemake with the denoise rule as the target, changing the number of cores to match your machine:
```
#snakemake --use-conda dada2_pe_denoise --cores 4
snakemake --use-conda dada2_pse_denoise --cores 4

#change cores when you like
```


#Pausing after the denoise step allows you to make changes before proceeding:

#Check the table summaries and representative sequence lengths to determine if DADA2 or Deblur parameters need to be modified. If so, you can rename or delete the output directories and then rerun the denoise rule.

#View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters "alpha_max_depth" and "core_sampling_depth" in config.yaml.

#Decide whether to filter your feature table and representative sequences by taxonomy or feature ID. After the taxonomy step, you can examine the taxonomy summary and bar plot to aid your decision. If you do filter your data, all output from that point on will go in a separate folder 

#Wso you can compare output with and without filtering.
#Unfiltered mode
#Continue the workflow without filtering (for now). If you are satisfied with your parameters and files, run the taxonomy rule (for unfiltered data):

```
snakemake --use-conda dada2_pe_taxonomy_unfiltered --cores 4
```

#Next, run the diversity rule (for unfiltered data):

```
snakemake --use-conda dada2_pe_diversity_unfiltered --cores 4
```

Finally, run the report rule (for unfiltered data):

```
snakemake --use-conda dada2_pe_report_unfiltered --cores 4

```





#Filtered mode
#After viewing the unfiltered results—the taxonomy summary and taxa barplot, the representative sequence summary plot and table, or the list of unassigned and potential outlier representative sequences—the user may wish to filter (remove) certain taxonomic groups or representative sequences. If so, the user should first check the following parameters and/or files:

#copy 2-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv to 00-data/repseqs_to_filter_dada2-pe.tsv to filter outliers, or manually include feature IDs in 00-data/repseqs_to_filter_dada2-pe.tsv to filter those feature IDs (change "dada2-pe" to "dada2-se" or "deblur-se" as appropriate);

#exclude_terms in config.yaml – add taxa to exclude from representative sequences, if desired;
repseq_min_length and repseq_max_length in config.yaml – set minimum and/or maximum lengths for filtering representative sequences, if desired;

#repseq_min_abundance and repseq_min_prevalence in config.yaml – set minimum abundance and/or prevalence values for filtering representative sequences, if desired.

#Now we are ready to filter the representative sequences and feature table, generate new summaries, and generate a new taxonomy bar plot, by running the taxonomy rule (for filtered data):

```
snakemake --use-conda dada2_pe_taxonomy_filtered --cores 4
```

Next, run the diversity rule (for filtered data):


```
snakemake --use-conda dada2_pe_diversity_filtered --cores 4

```

#Finally, run the report rule (for filtered data):

```
snakemake --use-conda dada2_pe_report_filtered --cores 1

```


```
View output
View report and output files
Open your HTML report (e.g., 03-reports/report_dada2-pe_unfiltered.html) in Chrome{target="_blank"} or Firefox{target="_blank"}. To view the linked files:

QZV (QIIME 2 visualization): click to download, then drag and drop in https://view.qiime2.org{target="_blank"}. Empress trees (e.g., rooted_tree.qzv) may take more than 10 minutes to load.
TSV (tab-separated values): click to download, then open in Microsoft Excel or Tabview (command line tool that comes with Tourmaline).
PDF (portable document format): click to open and view in new tab.
Downloaded files can be deleted after viewing because they are already stored in your Tourmaline directory.

More tips
Troubleshooting

#The whole workflow with test data should take ~3–5 minutes to complete. A normal dataset may take several hours to complete.

#If any of the above commands don't work, read the error messages carefully, try to figure out what went wrong, and attempt to fix the offending file. A common issue is the file paths in your FASTQ manifest file need to be updated.


#If you are running in a Docker container and you get an error like "Signals.SIGKILL: 9", you probably need to give Docker more memory. See the Wiki section on Installation options (https://github.com/aomlomics/tourmaline#:~:text=View%20output,the%20command%20below%3A).

Power tips

#The whole workflow can be run with just the command snakemake dada2_pe_report_unfiltered (without filtering representative sequences) or snakemake dada2_pe_report_filtered (after filtering representative sequences). Warning: If your parameters are not optimized, the results will be suboptimal (garbage in, garbage out).

#If you want to make a fresh run and not save the previous output, simply delete the output directories (e.g., 02-output-{method}-{filter} and 03-report) generated in the previous run. If you want to save these outputs and rerun with different parameters, you can change the name of the output directories and report files to something informative and leave them in the Tourmaline directory.

#You can always delete any file you want to regenerate. Then there are several ways to regenerate it: run snakemake FILE and Snakemake will determine which rules (commands) need to be run to generate that file; or, run snakemake RULE where the rule generates the desired file as output.

#If you've run Tourmaline on your dataset before, you can speed up the setup process and initialize a new Tourmaline directory (e.g., tourmaline-new) with the some of the files and symlinks of the existing one (e.g., tourmaline-existing) using the command below:
```
```
cd /path/to/tourmaline-new
scripts/initialize_dir_from_existing_tourmaline_dir.sh /path/to/tourmaline-existing
```


```
#send env to scratch
rsync -av /home/username/miniconda3/envs/myenv /scratch/username/myenv

```

```
clean conda env: 
conda clean --all
```

#tools
conda create -n bioinfo_env -c conda-forge -c bioconda multiqc fastqc samtools


multiqc: "/scratch/user/obayomi/.conda/envs/bioinfo/bin/multiqc"
    fastqc: "/scratch/user/obayomi/.conda/envs/bioinfo/bin/fastqc"
    parallel: "/scratch/user/obayomi/.conda/envs/bioinfo/bin/parallel"
    run_pear:  "pear" #"/scratch/user/obayomi/projects/qiime2/run_pear.pl"

myqueue

#https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/020-tutorial-upstream/020-metadata.html qiime2

snakamak conda env  for —sankemake-workflow-qiime2 pipeline



#https://github.com/SilasK/16S-dada2 : conda env dada2_env


https://github.com/a-h-b/dadasnake : dadasnake



```
#
https://www.youtube.com/embed/_wUGzqEjg6A

consider you ran a previous env
conda activate env
#export the env
conda env export > snakemake_amplicon.yml 
conda env create -f snakemake_amplicon.yml 
conda create --name snakemake --clone snakemake_amplicons_workflow
```

```
# Initialize Git (If Not Already Done

git init

#Connect to a GitHub Repository (Optional)
git remote add origin https://github.com/kaboyo/codes.git

#Verify the remote link:
git remote -v

#Stage, Commit, and Push Files
git add .

#Commit the Changes:

```


