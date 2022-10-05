[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/commit-activity)
[![GitHub contributors](https://img.shields.io/github/contributors/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/contributors/)
[![GitHub issues](https://img.shields.io/github/issues/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/issues/)

<p align="center"><img width="528" alt="MAW" src="https://user-images.githubusercontent.com/30716951/168855653-ae2efaa1-cbaf-4215-a04e-13bcd88ac46f.png"></p>


This repository hosts Metabolome Annotation Workflow (MAW). The workflow takes MS2 .mzML format data files as an input in R. It performs spectral database dereplication using R Package [Spectra](https://rformassspectrometry.github.io/Spectra/) and compound database dereplication using [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/). Final candidate selection is done in Python using [RDKit](https://www.rdkit.org/) and [PubChemPy](https://pubchempy.readthedocs.io/en/latest/).The classification of the tentative candidates from the input data are classified using [CANOPUS](https://bio.informatik.uni-jena.de/software/canopus/) and [ClassyFire](http://classyfire.wishartlab.com/).

## Install MAW

Install Docker on your [MAC OS](https://www.docker.com/get-started/) OR on [Linux](https://docs.docker.com/engine/install/ubuntu/). <br>

Since MAW is implemented in R and Python so we have two separate Docker images. 

### Pull maw-r Docker Image
To create a R-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-r
```
This will creat a R-Docker image on your system. This image contains /opt/workdir as the working directory which contains the following files and folders:
1. /data (folder)
2. /data/example_Tyrosine.mzML
3. /data/hmdb.rda
4. /data/gnps.rda
5. /data/mbankNIST.rda
6. Workflow_R_Functions.r (R function script)
7. Workflow_R_Script.r (R example script)
8. install_packages.R (R package installation script, these pacakges are already installed the docker container) <br>

### Pull maw-py Docker Image
To create a Python-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-py
```
This will creat a R-Docker image on your system. This image contains /opt/workdir as the working directory which contains the following files and folders:
1. /data (folder) 
2. /data/example_Tyrosine/ (folder)
3. /data/hmdb_dframe_str.csv
4. Workflow_Python_Functions.py (Python function script)
5. Workflow_Python_Script.py (Python example script)
6. requirements.txt (Python package installation script, these pacakges are already installed the docker container) <br>

## Use Case

1. First Run maw-r:

```
docker run --name example_maw-r -i -t zmahnoor/maw-r /bin/bash
```
Once you enter the docker container, you can run the script and check the results.
```
Rscript --no-save --no-restore --verbose Workflow_R_Script.r >outputFile.txt 2>&1
```
This command will run the Workflow_R_Script.r which is an example script for /data/example_Tyrosine.mzML.

2. Run maw-py
```
docker run --name example_maw-py -i -t zmahnoor/maw-py /bin/bash
```
Once you enter the docker container, you can run the script and check the results.
```
python3 Workflow_Python_Script.py
```
This command will run the Workflow_Python_Script.py which is an example script for /data/example_maw-r_results/. In order to leave the container without killing the process, add $ at the end of commands inside of the container. Then disown the PID and leave the container with CTRL+P and CTRL+Q.


For details on using the workflow on Jupyter notebooks in a more interactive mode, please follow the Tutorial part on wiki page of this repository

## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
