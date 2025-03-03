#!/bin/bash

# Example Snakefile for 16S Amplicon Sequencing Pipeline

# Define input/output files
FASTQ_DIR = "data/raw_fastq"
OUTPUT_DIR = "results"
TRIMMED_DIR = "results/trimmed"
MERGED_DIR = "results/merged"
OTU_DIR = "results/otu_table"

# List of samples (adjust according to your files)
SAMPLES = ["sample1", "sample2", "sample3"]

# Rule to run FastQC
rule fastqc:
    input:
        fq1 = f"{FASTQ_DIR}/{sample}_R1.fastq",
        fq2 = f"{FASTQ_DIR}/{sample}_R2.fastq"
    output:
        "results/fastqc/{sample}_R1_fastqc.html",
        "results/fastqc/{sample}_R2_fastqc.html"
    shell:
        "fastqc {input.fq1} {input.fq2} -o results/fastqc"

# Rule for trimming with Cutadapt
rule trim_reads:
    input:
        fq1 = f"{FASTQ_DIR}/{sample}_R1.fastq",
        fq2 = f"{FASTQ_DIR}/{sample}_R2.fastq"
    output:
        trimmed_fq1 = f"{TRIMMED_DIR}/{sample}_R1_trimmed.fastq",
        trimmed_fq2 = f"{TRIMMED_DIR}/{sample}_R2_trimmed.fastq"
    shell:
        "cutadapt -q 20 -o {output.trimmed_fq1} -p {output.trimmed_fq2} {input.fq1} {input.fq2}"

# Rule to merge paired-end reads (using VSEARCH)
rule merge_reads:
    input:
        fq1 = f"{TRIMMED_DIR}/{sample}_R1_trimmed.fastq",
        fq2 = f"{TRIMMED_DIR}/{sample}_R2_trimmed.fastq"
    output:
        merged = f"{MERGED_DIR}/{sample}_merged.fastq"
    shell:
        "vsearch --fastq_mergepairs {input.fq1} --reverse {input.fq2} --fastqout {output.merged}"

# Rule for OTU picking (using VSEARCH)
rule otu_picking:
    input:
        merged_fq = f"{MERGED_DIR}/{sample}_merged.fastq"
    output:
        otu_table = f"{OTU_DIR}/{sample}_otu_table.txt"
    shell:
        "vsearch --cluster_size {input.merged_fq} --id 0.97 --centroids {output.otu_table}"

# Optional: Rule for taxonomic classification (using QIIME2 or VSEARCH)
rule classify_taxonomy:
    input:
        otu_table = f"{OTU_DIR}/{sample}_otu_table.txt"
    output:
        taxonomic_assignment = f"{OTU_DIR}/{sample}_taxonomy.txt"
    shell:
        "vsearch --usearch_global {input.otu_table} --db reference_database.fasta --id 0.97 --blast6out {output.taxonomic_assignment}"

# Final rule: combine results or additional analysis steps
rule combine_results:
    input:
        expand("{OTU_DIR}/{sample}_otu_table.txt", sample=SAMPLES)
    output:
        "results/combined_otu_table.txt"
    shell:
        "cat {input} > {output}"

