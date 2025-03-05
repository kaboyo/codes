#brief guide https://www.youtube.com/embed/_wUGzqEjg6A
import os
import glob

SRA,FRR = globe_wildcards (rawReads/(sra)_[frr].fastq.gz)

rule all:
    input: 
        expand("rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA), frr=FRR, extension=["html", "zip"])

     # TODO: add your rules for other steps here
rule rawFastqc:
    input:  
        rawReads="rawReads/(sra)_[frr].fastq.gz"
        output:
            zip="rawFastqc/ {sra}_{ffr}_fastqc.zip",
            html="rawFastqc/{sra}_{frr}_fastqc_report.html"
            threads:
                2
                params:
                    path="rawFastqc/"
                    
                shell:
                """
                fastqc{input.rawread} --threads {threads}  -o {params.path}
                """

    rule trimmomaticd:
    input:
        read1= "rawReads/(ids)_1.fastq.gz",
        read2= "rawReads/(ids)_2.fastq.gz"
        
        output: 
            forwardpaired="trimmedReads/{ids}_1P.fastq",
            reversepaired="trimmedReads/{ids}_2P.fastq",

            threads: 
                4
                params:
                    basename="trimmedReads/{sra}.fastq",
                    log="trimmedReads/{sra}.log"
                shell: 
                """
                {output.trimmed_fq1} {output.unpaired_fq1} {output.trimmed_fq2} {output.unpaired_fq2} ILLUMINACLIP:/path/to/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
                """
    rule cutadapter:
        input:
            read1= "rawReads/(ids)_1.fastq.gz",
            read2= "rawReads/(ids)_2.fastq.gz"
            output:
                forwardpaired="trimmedReads/{ids}_1P.fastq",
                reversepaired="trimmedReads/{ids}_2P.fastq",
                threads: 
                4
                params:
                    basename="trimmedReads/{sra}.fastq",
                    log="trimmedReads/{sra}.log"
                shell: 
                """
                cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {output.forwardpaired} -p {output.reversepaired} {input.read1} {input.read2}
                ""