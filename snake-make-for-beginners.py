#http://ivory.idyll.org/blog/2023-snakemake-slithering-section-1.html
#snakemake for doing bioinformatics - a beginner's guide (part 1)
#Installation and setup!

```
#Setup and installation

#I suggest working in a new directory.
#You'll need to install snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) & 
# sourmash (https://sourmash.readthedocs.io/en/latest/#installing-sourmash). 
# We suggest using mamba, via miniforge/mambaforge (https://github.com/conda-forge/miniforge#mambaforge), for this.



```
scp eo2r24@iridis5.soton.ac.uk:/scratch/eo2r24/SNAKEMAKE/tourmaline/00-data/metadata.tsv /Users/expeditoolimi/Documents/PRJNA358488/codes 
```



#https://youtu.be/ow0nZsbgusw?t=316
#https://youtu.be/SMDdfcpONX8?t=64 snakemake 

#Snakemake rule all

#Hi everyone, and welcome back to the Snakemake series
#In this video, I will show you how Snakemake decides what rules to run and when, and also why we need a rule all that has no output and only inputs.
#Let's start with a simple rule that copies file a and renames it file b. If we launch the pipeline, we have only the copyab rule run, 
#which generates the b file. Now let's add a rule that makes a copy of file b and renames it file c. 

rule copyab:
  input: "a.txt"
  output: "b.txt"
  shell: "cp {input} {output}"

rule copybc:
  input: "b.txt"
  output: "{x}.txt"
  shell: "cp {input} {output}"

#If we launch the pipeline, we would expect to get file c, but this does not happen. Only the copyab rule is run instead. 
# This is because Snakemake takes the first rule and tries to get its output, in this case, file b. Since to get file b, we just need 
# file a and file a is already present, the only rule that needs to be run is rule copyab. Any other rule below is irrelevant.

#Let's see what happens if we move the rule copybc to the top. Now when we launch the pipeline, both rules are run, and the file c is generated. 
# What is going on is that Snakemake tries to run the first rule, but the b input file is absent, so it looks for another rule that can generate 
# the file and finds rule copyab that generates the b file and runs it. 

rule copybc:
  input: "b.txt"
  output: "{x}.txt"
  shell: "cp {input} {output}"

rule copyab:
  input: "a.txt"
  output: "b.txt"
  shell: "cp {input} {output}"

#After running copyab, the copybc is launched, and we end up running both rules. So far, so good, if we place the rules in the correct order, 
# we can get them all to run. But still, we can run into problems with divergent branches. For example, if we want to generate a second copy of 
# file b and name it d, we need to write a new rule. However, since the first rule of the pipeline does not depend on the new rule, 
# the rule is not run. 

rule copybc:
  input: "b.txt"
  output: "c.txt"
  shell: "cp {input} {output}"

rule copybd:
   input: "b.txt"
   output: "d.txt"
   shell: "cp {input} {output}"

rule copyab:
  input: "a.txt"
  output: "b.txt"
  shell: "cp {input} {output}"

# One way to solve this problem would be to come up with a rule that needs both files, c, and d, as inputs. 
# For example, we can write the following rule catcd. Now if we launch the pipeline, both files c and d are generated together with this 
# additional file x.

rule catcd:
  input: "c.txt", "d.txt"
  output:  "x.txt"
  shell: "cat {input} x {output}"

#But in this case, I just want files c and d. I do not care about x. Let's see what happens if we only use the input for rule catcd. 
# Oddly enough, the pipeline launches and generates the files as expected. If we check the run, all rules are run normally, 
# but the catcd one is flagged as local rule, and no actual shell command is associated with it. So for rule catcd, we can leave only the input, 
# and we are good. 

rule catcd:
  input: "c.txt", "d.txt"
output:  "x.txt"
shell: "cat {input} x {output}"

#We will rename also the rule, rule all because that is the general convention. 
# And this is the reason why the output of the whole pipeline is specified as input within the rule all and not as output, 
# which initially can result counterintuitive. Now we can give any order to the other rules since this does not longer matter. 
# My preference is to write them in the order they are executed. 

rule all:
  input: "c.txt", "d.txt"

#Well, that's all you needed to know to understand rule all and how Snakemake decides what rules to run.


#Building Bioinformatics Pipelines with Snakemake
https://mtbgenomicsworkshop.readthedocs.io/en/latest/material/day4/reproducible_research.html