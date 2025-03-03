import os
import glob

SRA = glob_wildcards("rawReads/{sra}.fastq.gz")

rule all:
    input: 
        expand("rawQC/{sra}_fastqc.{extension}", sra=SRA, extension=["html", "zip"])

rule rawFastqc:
    input:  
        rawReads="rawReads/{sra}.fastq.gz"
    output:
        zip="rawFastqc/{sra}_fastqc.zip",
        html="rawFastqc/{sra}_fastqc_report.html"
    threads: 2
    params:
        path="rawFastqc/"
    shell:
        """
        fastqc {input.rawReads} --threads {threads} -o {params.path}
        """

rule trimmomatic:
    input:
        read="rawReads/{sra}.fastq.gz"
    output: 
        trimmed="trimmedReads/{sra}_trimmed.fastq.gz"
    threads: 4
    params:
        log="trimmedReads/{sra}.log"
    shell: 
        """
        trimmomatic SE -threads {threads} {input.read} {output.trimmed} \
        ILLUMINACLIP:/path/to/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > {params.log} 2>&1
        """

rule cutadapter:
    input:
        read="rawReads/{sra}.fastq.gz"
    output:
        trimmed="trimmedReads/{sra}_cutadapt.fastq.gz"
    threads: 4
    params:
        log="trimmedReads/{sra}_cutadapt.log"
    shell: 
        """
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -o {output.trimmed} {input.read} > {params.log} 2>&1
        """
