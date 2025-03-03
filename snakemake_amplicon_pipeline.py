import os
import glob

SRA, FRR = glob_wildcards("rawReads/{sra}_{frr}.fastq.gz")

rule all:
    input: 
        expand("rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR, extension=["html", "zip"])

rule rawFastqc:
    input:  
        rawReads="rawReads/{sra}_{frr}.fastq.gz"
    output:
        zip="rawFastqc/{sra}_{frr}_fastqc.zip",
        html="rawFastqc/{sra}_{frr}_fastqc_report.html"
    threads: 2
    params:
        path="rawFastqc/"
    shell:
        """
        fastqc {input.rawReads} --threads {threads} -o {params.path}
        """

rule trimmomaticd:
    input:
        read1="rawReads/{ids}_1.fastq.gz",
        read2="rawReads/{ids}_2.fastq.gz"
    output: 
        forwardpaired="trimmedReads/{ids}_1P.fastq",
        reversepaired="trimmedReads/{ids}_2P.fastq"
    threads: 4
    params:
        basename="trimmedReads/{ids}.fastq",
        log="trimmedReads/{ids}.log"
    shell: 
        """
        trimmomatic PE -threads {threads} {input.read1} {input.read2} \
        {output.forwardpaired} {output.reversepaired} \
        ILLUMINACLIP:/path/to/adapters/TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > {params.log} 2>&1
        """

rule cutadapter:
    input:
        read1="rawReads/{ids}_1.fastq.gz",
        read2="rawReads/{ids}_2.fastq.gz"
    output:
        forwardpaired="trimmedReads/{ids}_1P.fastq",
        reversepaired="trimmedReads/{ids}_2P.fastq"
    threads: 4
    params:
        basename="trimmedReads/{ids}.fastq",
        log="trimmedReads/{ids}.log"
    shell: 
        """
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o {output.forwardpaired} -p {output.reversepaired} {input.read1} {input.read2} \
        > {params.log} 2>&1
        """
