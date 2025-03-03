

#nextflow pipeline


```
#https://github.com/nf-core/ampliseq/blob/master/CITATIONS.md /details of the different pipeline packages

#Install pipeline

Check nfcore: https://nf-co.re/docs/usage/getting_started/installation

/1 Quick-start installation
/2 Bioconda installation [considered bioconda installation]

```

```
#To run a pipeline:

#All software dependencies must be installed (Nextflow + Docker / Singularity / Conda). See Dependency installation for more information.

#Run Nextflow’s “hello world” example to confirm all tools are installed and operational:

nextflow run hello

#Choose a pipeline to run. See the available nf-core pipelines. If you have nf-core tools installed, run nf-core pipelines list.

check nf-core/tools for for intallation of nf-core tools in the hood of the nextflow env (assuming that you installed nextflow as a conda environment (https://nf-co.re/tools)

#Configure Nextflow to run on your system:

1/The simplest way to run is with -profile docker (or singularity). This instructs Nextflow to execute jobs locally, with Docker (or Singularity) to fulfill software dependencies.

2/Conda is also supported with -profile conda. However, this option is not recommended, as reproducibility of results can’t be guaranteed without containerization.

3/If you are a member of one of the listed institutions, use the institutional config file created for your institution.

4/See Nextflow configuration for advanced Nextflow configuration options.

#Run the tests for your pipeline in the terminal to confirm everything is working:

nextflow run nf-core/<pipeline_name> -profile test,docker --outdir <OUTDIR>
Replace <pipeline_name> with the name of an nf-core pipeline.

#If you don’t have Docker installed, replace docker with either singularity or conda.

Nextflow will pull the code from the GitHub repository and fetch the software requirements automatically, so there is no need to download anything first.

#If the pipeline fails, see Troubleshooting or ask for help on the nf-core Slack channel for your pipeline.

#Read the pipeline documentation to see which command-line parameters are required. This will be specific to your data type and usage.


#download the nf-core pipeline that you are interested in running and have it locally available on the cluster as running on an hpc doesnot let you download staff 


```

#single cell
https://github.com/csgenetics/csgenetics_scrnaseq?tab=readme-ov-file#cs-genetics-scrna-seq-pipeline

https://github.com/Goldrathlab/Spatial-TRM-paper

https://github.com/maximilian-heeg/vizgen-segmentation/

kallisto:https://github.com/cbcrg/kallisto-nf-reproduce/tree/nbt-v1.0?tab=readme-ov-file

cancer genomics:https://github.com/abhinandan0y/CancerDiagnosisAI

https://github.com/selcukorkmaz/fastml


```
#Pull the nfcore-pipeline for amplicon sequencing
-activate nextflow and nf-core tools installed in the same hood
-call 
#nf-core pipelines list
-all available tools will be listed
-check the ampliseq tool
#-use the: nf-core download ampliseq; 
and follow the process to get singularity installed
-the program will be pulled and ready for running 

# command to pull the program
nf-core pipelines download ampliseq -r 2.12.0 -s singularity


#Viola!

#command [for single end]
NXF_VER=24.10.4 nextflow run nf-core/ampliseq -r 2.12.0 -profile conda --input "samplesheet.tsv" --skip_cutadapt --outdir amplis_test --single_end
#copy files;
scp -r /Users/expeditoolimi/Documents/NEXTFLOW/PRJNA358488 eo2r24@iridis5.soton.ac.uk:/scratch/eo2r24/nexflow/PRJNA358488 

scp -r eo2r24@iridis5.soton.ac.uk:/scratch/eo2r24/nexflow/PRJNA358488/ampli /Users/expeditoolimi/Documents/NEXTFLOW/PRJNA358488


#unzip the .gz
tar -xzvf filename.tar.gz

```


```
#Test pipeline using the conda profile but the locall installed packages
#worked: NXF_VER=24.10.4 nextflow run nf-core-ampliseq_2.12.0/2_12_0/main.nf -profile test,conda --outdir Test


#Failed
NXF_VER=24.10.4 nextflow run nf-core-ampliseq_2.12.0/2_12_0/main.nf -profile conda --input "samplesheet.tsv"  --single_end --outdir verena_test


```



https://oldsite.nf-co.re/launch?id=1740251247_05e28e4d951c
https://nf-co.re/ampliseq/2.12.0/docs/usage/
https://slurm.schedmd.com/quickstart.html
https://nf-co.re/docs/nf-core-tools/installation
https://nf-co.re/docs/usage/getting_started/configuration

