---
layout: post
title: Python on HPC3
subheading: How to use python/R/other programming language on the HPC3
author: Emmanuel Dollinger
categories: hpc
banner: /assets/images/other/anaconda.jpeg
tags: hpc python
sidebar: []
usemathjax: true
---

Author: Emmanuel Dollinger

## How to run python on the HPC3

Sam Morabito recently showed me a trick to get python working on the HPC3. This makes running python on the HPC3 actually useable. This trick relies on using conda environments, you can read up on conda [here](https://conda.io/projects/conda/en/latest/user-guide/index.html) and on conda virtual environments [here](https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html).

## Scope:

This will only work if you are running your code in a programming language such as python or R. If you are running a package directly from bash such as cellranger this isn't what you want to do.

I won't go over basic HPC3/bash code here, see [this link](https://rcic.uci.edu/hpc3/slurm.html) for HPC3 stuff.

We will install Scanpy, which is a nice single cell sequencing package for python. This is a nice usecase because Scanpy is both very common and not in conda.

## Create new conda environment

```bash

ssh -y edolling@hpc3.rcic.uci.edu

module load anaconda

conda init #do whatever it tells you to do, probably will have to quit and login in again to the HPC3.

conda create -n EnvName python=3.6 # create a new environment with specified python version

```

It will ask you to install packages, continue through.

Now we have a virtual conda environment, that we can install packages in directly. This totally circumvents the onerous requirement to have modulefiles and such for each package. Best practices are to create a virtual environment for each project.

## Install packages that are in conda

Every time you want to use a virtual environment, you need to activate it (again see the env tutorial above).

```bash

conda activate EnvName

```

The environment name should change from (base) to (EnvName).

Best practice for packages that are not in conda is to first install everything that is in conda and then download via pip. The reason for this is conda can update its packages but not other packages, and doesn't "see" pip updating packages. See [this post](https://www.anaconda.com/blog/using-pip-in-a-conda-environment) for more.

```bash

conda install pandas

# Install a bunch of packages

conda install matplotlib

# Install moar packages

conda install seaborn

# INSTALL ALL THE PACKAGES

# etc

```

If you already have a conda env that has all the packages you want, you can do:

```bash

conda activate OldEnv

conda list -e > packagelist.txt

conda create -n NewEnv --file packagelist.txt

```

You can also create the conda env and pass the file to `conda install` afterwards. NB this will error if you have pip installed packages in the old conda env.

## Install packages not in conda

Once you've installed most/all of the packages in conda, you can simply pip install scanpy.

```bash

conda activate EnvName # Don't forget to do this

pip install scanpy

```

## Load scanpy in python file

Now we will submit a job to the scheduler that just loads scanpy and quits. You need two scripts, a slurm script that submits the python script and the python script.

slurm.sub script:

```bash
#!/bin/bash

#SBATCH --job-name=scanpytest      ## Name of the job.
#SBATCH -A qnie_lab     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=../out ## out log file

# Run the following two lines every time you submit a python script to slurm, this tells slurm about your conda env and loads it.
source ~/.bashrc

conda activate trVAE

# This next line just runs the python script

python3 testscipt.py

```

testscipt.py script:

```python
import scanpy as sc

print("Scanpy successfully loaded.")

```

In ../out:
```bash
Scanpy successfully loaded.

```

That's it! NB this will also work with R and any other language that conda supports.
