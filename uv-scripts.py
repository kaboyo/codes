#prepared, thanks to Stephen Turner blogs
#https://docs.astral.sh/uv/
#https://blog.stephenturner.us/p/uv-part-1-running-scripts-and-tools?r=wxa1&utm_campaign=post&utm_medium=web
#https://blog.stephenturner.us/p/uv-part-1-running-scripts-and-tools?r=wxa1&utm_campaign=post&utm_medium=web
#uv, part 1: running scripts and tools

#Running scripts
#One of the most immediately useful things you can do with uv is running Python scripts. 
# With no dependencies you can do this with: uv run example.py. 
# This doesn’t do anything more than running python example.py would do.
```
print("Hello world")
```

```
uv run example.py
Hello world
```


#But what if we have imports that aren’t part of the standard library? 
# Here’s a script that uses Biopython to reverse complement a DNA sequence, 
# and prints the original and the reverse complement.

```
from Bio.Seq import Seq

my_seq = Seq("GATTACA")
print(my_seq)

rc = my_seq.reverse_complement()
print(rc)
```

#If we just try to run that with python example.py or uv run example.py, 
#it will fail because we haven’t installed Biopython.

```
$ uv run example.py
Traceback (most recent call last):
  File "/Users/turner/Downloads/uv/example.py", line 1, in <module>
    from Bio.Seq import Seq
ModuleNotFoundError: No module named 'Bio'
```

#You could go through the trouble of setting up a new Conda environment or venv, 
# the pip installing biopython, but this adds extra maintenance burden and time solving your Conda environment.

#Running scripts with inline metadata
#PEP 723 specifies a metadata format that can be embedded in 
# single-file Python scripts to declare dependencies, and uv 
# supports declaring dependencies inline according to PEP 723.

#Here’s the same script with the dependencies declared inline:
uv add --script example.py biopython

```

#This modifies the script in place, adding the biopython dependency in this specialized form of a comment at the top of the script.

#insert script here
```
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython",
# ]
# ///
from Bio.Seq import Seq

my_seq = Seq("AGTACACTGGT")
print(my_seq)

rc = my_seq.reverse_complement()
print(rc)
```

#We can now run it with uv run, and when we do so, uv understands this inline metadata at 
# the top and creates this temporary virtual environment on the fly, installs Biopython, 
# and runs the code. It does this all in less than one second (package installation took only 10 milliseconds) 🔥.

```
$ time uv run example.py
Reading inline script metadata from `example.py`
Installed 2 packages in 10ms

GATTACA
TGTAATC

real	0m0.767s
user	0m0.223s
sys	0m0.231s
```


#That on-demand environment is cached, so running it again takes milliseconds 🔥.

```
time uv run example.py
Reading inline script metadata from `example.py`

GATTACA
TGTAATC

real	0m0.080s
user	0m0.048s
sys	0m0.027s
```

#Executable scripts with a uv shebang

#You can add a shebang line to the top of your script that’ll allow you to run it as an executable:
#  #!/usr/bin/env -S uv run. 
# I also added the --quiet flag so you don’t get any additional output (“Reading inline script metadata…”). 
# I learned about this on Simon Willison’s blog. 


# Here’s the script with the shebang:
```
#!/usr/bin/env -S uv run --quiet

# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython",
# ]
# ///

from Bio.Seq import Seq

my_seq = Seq("AGTACACTGGT")
print(my_seq)

rc = my_seq.reverse_complement()
print(rc)
```

#You can chmod 755 your script and run it as an executable on any machine where uv is available, 
# and if the dependencies you specify with inline metadata will automatically be installed 
# and cached into a temporary isolated environment, using the correct version of Python as specified.




#Running tools
#From PyPI
#uv can also run tools from PyPI. Here’s an example using the tool seqkit:

#With uv you can also use tools that are part of Python packages. 
# You can do this with uv tool run <toolname> or the shortcut, 
# uvx <toolname> (uvx is meant to be a faster drop-in replacement for pipx). 
# At my previous job we wrote a little Python package called vcferr (paper, PyPI, GitHub) 
# that takes a VCF and changes some of the genotypes in ways you control. 
# You could set up a virtual environment and pip install it, or a conda environment and 
# conda install it. But you can avoid creating and managing these environments and 
# just use uv to run the tool directly. Because vcferr is on PyPI, we can just use uv tool run vcferr, 
# or shorter, uvx vcferr. You can see that installation of its dependencies into a temporary environment 
# took only 5 milliseconds, and the entire process required only 3 seconds.
# Because this is now cached, running it again would only take milliseconds 🔥.

```
time uvx vcferr --help
Installed 3 packages in 5ms

Usage: vcferr [OPTIONS] <input_vcf>

Options:
  -s, --sample TEXT        ID of sample in VCF file 
                           to be simulated[required]
  -o, --output_vcf TEXT    Output VCF file containing simulated 
                           genotypes ex:example.sim.vcf.gz
  ...(truncated)...
  -a, --seed INTEGER       Random number seed
  --help                   Show this message and exit.

real	0m3.351s
user	0m0.280s
sys	0m0.184s
```

#You could use this for little tools where you don’t want to set up a full environment, such as running a 
# code formatter/linter like black (uvx black) or ruff (uvx ruff). If you use the tool a lot, you could use 
# uv tool install <toolname> to make the tool available in a special directory on your PATH.

#From GitHub

#We can also run a tool directly from GitHub. A few months ago 
# I wrote about making a Python CLI using click, and making this a package using a cookiecutter template.
#link here:https://blog.stephenturner.us/p/python-cli-click-cookiecutter


#This tool, caffeinated (GitHub, PyPI), is a silly little command line tool to tell you how much caffeine 
# remains in your system at bedtime based on how much you consume in the morning. It’s on PyPI so you 
# could just use uvx caffeinated to run it. However, what if it wasn’t on PyPI? We can run the tool directly 
# from GitHub by passing in the --from option and giving it the GitHub URL. You can see that it installed 
# the two required dependencies in 1 millisecond, and runs the code in 1.6 seconds. The tool is cached so 
# running it a second time takes only milliseconds 🔥.

```
$ time uvx --from \
  git+https://github.com/stephenturner/caffeinated \
  caffeinated -c 200 -s 0600 -b 2000

 Updated https://github.com/stephenturner/caffeinated (762d601)
   Built caffeinated @ git+https://github.com/stephenturner/caffeinated@762d60111ac84f4af5031141ea29232d26659f9f
Installed 2 packages in 1ms

You would have 39.7mg of caffeine in your system 
if you went to bed at 8:00pm (in 14.0 hours).
That's like having 44% of a cup of coffee before bed.

real	0m1.678s
user	0m0.611s
sys	0m1.481s
```


#Tools with many dependencies
#In the examples above the tools had a very lightweight dependency stack. What about a tool with more first and 
# second order dependencies? The nf-core/tools Python package has a few dependencies, each of those has a few 
# second-order and deeper dependencies. When we run it the first time with uvx nf-core it takes ~7.5 seconds to 
# resolve the dependency stack of 95 packages, but less than one second to install those dependencies 🔥.

```
$ time uvx nf-core
Installed 95 packages in 255ms
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.1.1 - https://nf-co.re

 Usage: nf-core [OPTIONS] COMMAND [ARGS]...

(...truncated...)

real	0m7.501s
user	0m2.089s
sys	0m5.420s
```

#Since the result is cached, running it a second time required less than one second 🔥.


#Tools requiring other packages in the environment
#There’s a tool on PyPI called voila (documentation) that turns Jupyter notebooks into standalone web applications.
# Let’s take a look at this notebook that uses bqplot to create an interactive graphic. 
# Here’s what it looks like in Jupyter lab (I’ll cover running Jupyter with uv in another post).
#We can use voila to turn that notebook into an interactive web app. You may think of running something 
# like uvx voila bqplot.ipynb. But this won’t work — voila runs without an issue but it’ll complain that you 
# need a jupyter kernel, numpy, and bqplot installed to actually run the conversion. Because none of these are 
# dependencies of voila itself, uv happily runs voila. But because the thing you’re doing with voila requires 
# additional packages available that aren’t dependencies of voila itself, you’ll have to add the --with option, 
# specifying which additional packages you want installed in the on-demand environment you create when using uvx voila.

```
#We can use voila to turn that notebook into an interactive web app. You may think of running something
uvx --with jupyter,numpy,bqplot voila bqplot.ipynb
```
#The result is below. The first run with a clean cache installs 105 packages from PyPI, and took ~6 seconds.
# Subsequent runs using uv’s cache initialized in <1 second 🔥.


#More info
#The best way to learn more is by going through the Getting started section and the Guides sections on 
# the official docs at docs.astral.sh/uv.

#Arjan’s video on uv provides a good overview of what you can do with uv. https://youtu.be/qh98qOND6MI 


#Appreciation to Stephen Turner for the blog post.




#uv manages project dependencies and environments, with support for lockfiles, workspaces, and more, similar to rye or poetry:

```
uv init example-uv


cd example

uv add ruff

#Creating virtual environment at: .venv
#Resolved 2 packages in 374ms
#Prepared 1 package in 691ms
#Installed 1 package in 41ms
# + ruff==0.9.9



uv run ruff check
#All checks passed!

uv lock
#Resolved 2 packages in 2ms

uv sync
#Resolved 2 packages in 1ms
#Audited 1 package in 0.47ms
```

#Scripts
#uv manages dependencies and environments for single-file scripts.
#uv can run scripts with dependencies declared inline, similar to PEP 723:
uv init SCRIPTS
 cd SCRIPTS/
 echo 'import requests; print(requests.get("https://astral.sh"))' > example.py
 nano example.py 
 uv add --script example.py requests

```
#Then, run the script in an isolated virtual environment:
#detailed gude to scripts with uv and uvx: https://docs.astral.sh/uv/guides/scripts/#running-a-script-without-dependencies
```
uv run example.py

```

#Tools
#uv executes and installs command-line tools provided by Python packages, similar to pipx.
#Run a tool in an ephemeral environment using uvx (an alias for uv tool run):

```
uvx pycowsay 'hello world!'
#hello world!
```
Install a tool with uv tool install:#https://docs.astral.sh/uv/guides/tools/#installing-a-tool
uv tool install ruff
#Installed ruff==0.9.9
```

#Python versions
#uv installs Python and allows quickly switching between versions.
#Install multiple Python versions:
```
uv python install 3.10 3.11 3.12
#Installed 2 versions in 5.04s
# + cpython-3.10.16-macos-x86_64-none
# + cpython-3.12.9-macos-x86_64-none
```
#Download Python versions as needed:
```
uv venv --python 3.12.0
#Using CPython 3.12.0
#Creating virtual environment at: .venv
#Activate with: source .venv/bin/activate
uv run --python pypy@3.8 -- python
```
#Use a specific Python version in the current directory:

```
uv python pin 3.11
#Pinned `.python-version` to `3.11`
```

#The pip interface
#uv provides a drop-in replacement for common pip, pip-tools, and virtualenv commands.
#uv extends their interfaces with advanced features, such as dependency version overrides, platform-independent resolutions, reproducible resolutions, alternative resolution strategies, and more.
#Migrate to uv without changing your existing workflows — and experience a 10-100x speedup — with the uv pip interface.

#Compile requirements into a platform-independent requirements file:

```
uv pip compile docs/requirements.in  --universal --output-file docs/requirements.txt

or
uv pip compile docs/requirements.in \
   --universal \
   --output-file docs/requirements.txt
#Create a virtual environment:
uv venv
#Install the locked requirements:

uv pip sync docs/requirements.txt

```

#https://docs.astral.sh/uv/ for details on uv

#Creating a new project
#https://docs.astral.sh/uv/guides/projects/#python-version
#You can create a new Python project using the uv init command:
```
uv init hello-world
 cd hello-world/
 uv init
```

```
#Managing dependencies
uv add requests

#To remove a package, you can use uv remove:
uv remove requests
#To upgrade a package, run uv lock with the --upgrade-package flag:
uv lock --upgrade-package requests
#To install the locked dependencies, run uv sync:

```


#Running commands
#uv run can be used to run arbitrary scripts or commands in your project environment.
#Prior to every uv run invocation, uv will verify that the lockfile is up-to-date with the pyproject.toml, 
# and that the environment is up-to-date with the lockfile, keeping your project in-sync without the need for manual intervention. 
# uv run guarantees that your command is run in a consistent, locked environment.

#For example, to use flask:

```
uv add flask
uv run -- flask run -p 3000
```
#https://docs.astral.sh/uv/guides/projects/#next-steps

#to check> Python for R users:https://blog.stephenturner.us/p/python-for-r-users 
to check on Click: https://click.palletsprojects.com/en/stable/
to check on Cookiecutter: https://cookiecutter.readthedocs.io/en/stable/
panExplorer_workflow: https://github.com/SouthGreenPlatform/PanExplorer_workflow 

#install uv with pip
pip install uv
