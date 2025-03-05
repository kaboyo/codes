import os
import snakemake

# Define input data directory
RAW_READS = "data/raw_reads"
TRIMMED_READS = "data/trimmed_reads"
DENOISED_READS = "data/denoised_reads"
TAXONOMY = "data/taxonomy"

rule all:
    input:
        expand("results/otu_table_{{sample}}.qza", sample=glob_wildcards(os.path.join(RAW_READS, "SRR*.fq.gz")).sample),
        expand("results/taxonomy_{{sample}}.tsv", sample=glob_wildcards(os.path.join(RAW_READS, "SRR*.fq.gz")).sample)

rule fastqc:
    input:
        raw_reads = os.path.join(RAW_READS, "{sample}.fq.gz")
    output:
        qc_report = os.path.join("results/qc", "{sample}_fastqc.html")
    shell:
        "fastqc {input.raw_reads} -o results/qc"

rule trim_reads:
    input:
        raw_reads = os.path.join(RAW_READS, "{sample}.fq.gz")
    output:
        trimmed_reads = os.path.join(TRIMMED_READS, "{sample}_trimmed.fq.gz")
    shell:
        "trimmomatic SE {input.raw_reads} {output.trimmed_reads}"

rule dada2_denoise:
    input:
        trimmed_reads = os.path.join(TRIMMED_READS, "{sample}_trimmed.fq.gz")
    output:
        denoised_table = os.path.join(DENOISED_READS, "{sample}_table.qza")
    shell:
        "qiime dada2 denoise-single --i-demultiplexed-seqs {input.trimmed_reads} --o-table {output.denoised_table}"

rule assign_taxonomy:
    input:
        denoised_table = os.path.join(DENOISED_READS, "{sample}_table.qza")
    output:
        taxonomy = os.path.join(TAXONOMY, "{sample}_taxonomy.qza"),
        taxonomy_tsv = os.path.join("results", "taxonomy_{sample}.tsv")
    shell:
        "qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads {input.denoised_table} --o-classification {output.taxonomy}"

rule export_taxonomy:
    input:
        taxonomy = os.path.join(TAXONOMY, "{sample}_taxonomy.qza")
    output:
        taxonomy_tsv = "results/taxonomy_{sample}.tsv"
    shell:
        "qiime tools export --input-path {input.taxonomy} --output-path results && mv results/taxonomy.tsv {output.taxonomy_tsv}"



# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "glob",
#     "glob-wildcards",
#     "globe",
#     "os",
#     "snakemake",
# ]
# ///

import os
from snakemake.io import glob_wildcards

# Define input data directory
RAW_READS = "data/raw_reads"

# Automatically find all sample names that match the pattern
SAMPLES, = glob_wildcards(os.path.join(RAW_READS, "SRR55274*_R1.fq.gz"))

rule all:
    input:
        expand("results/qc/{sample}_fastqc.html", sample=SAMPLES)

rule fastqc:
    input:
        raw_reads = os.path.join(RAW_READS, "{sample}_R1.fq.gz")
    output:
        qc_report = os.path.join("results/qc", "{sample}_fastqc.html")
    shell:
        "fastqc {input.raw_reads} -o results/qc"
        
        
        
        
        
        
        
        
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "glob",
#     "glob-wildcards",
#     "os",
#     "snakemake",
# ]
# ///

import os
from snakemake.io import glob_wildcards

# Define input data directory
RAW_READS = "data/raw_reads"

# Use glob_wildcards to extract sample names from the raw reads
SAMPLES, = glob_wildcards(os.path.join(RAW_READS, "SRR55274*_R1.fq.gz"))

rule all:
    input:
        expand("results/qc/{sample}_fastqc.html", sample=SAMPLES)

rule fastqc:
    input:
        raw_reads = os.path.join(RAW_READS, "{sample}_R1.fq.gz")
    output:
        qc_report = os.path.join("results/qc", "{sample}_fastqc.html")
    shell:
        "fastqc {input.raw_reads} -o results/qc"

#create a sample table like this one. You can use the script prepare_sample_table.py for it. 
# The scripts searches for fastq(.gz) files inside a folder (structure). If you have paired end files they should have R1/R2 somewhere in the filename. If might be a good idea to simplify sample names.

```

./prepare_sample_table.py path/to/fastq(.gz)files

```