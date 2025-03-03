#How to create an env from the previous one (cloning)
#consider you ran a previous env
```
conda activate env
```
#export the env
```
conda env export > snakemake_amplicon.yml 
conda env create -f snakemake_amplicon.yml 
```
or clone  an env and rename it
```
conda create --name snakemake --clone snakemake_amplicons_workflow
```

#the code for snakemake_amplicons
Link: https://github.com/a-h-b/dadasnake 

#Tutorial
dadasnake is a Snakemake workflow to process amplicon sequencing data, from raw fastq-files to taxonomically assigned "OTU" tables, based on the DADA2 method. Running dadasnake could not be easier: it is called by a single command from the command line. With a human-readable configuration file and a simple sample table, its steps are adjustable to a wide array of input data and requirements. It is designed to run on a computing cluster using a single conda environment in multiple jobs triggered by Snakemake. dadasnake reports on intermediary steps and statistics in intuitive figures and tables. Final data output formats include biom format, phyloseq objects, and flexible text files or R data sets for easy integration in microbial ecology analysis scripts.


```
#Installing dadasnake

#For dadasnake to work, you need conda.
#Clone this repository to your disk:
git clone https://github.com/a-h-b/dadasnake.git

#Change into the dadasnake directory:
cd dadasnake
```

#if you want to submit the process running snakemake to the cluster:
```
cp auxiliary_files/dadasnake_allSubmit dadasnake
chmod 755 dadasnake
```

#if you want to keep the process running snakemake on the frontend using tmux:

```
cp auxiliary_files/dadasnake_tmux dadasnake
chmod 755 dadasnake
```
#Alternatively, if the above does not work, you can install a fixed snakemake version without mamba like so:
```
conda env create -f workflow/envs/snakemake_env.yml 
```

#initialize conda environments: This run sets up the conda environments that will be usable by all users
```
./dadasnake -i config/config.init.yaml 
```


#Optional test run: The test run does not need any databases. You should be able to start it by running

```
./dadasnake -l -n "TESTRUN" -r config/config.test.yaml
```

