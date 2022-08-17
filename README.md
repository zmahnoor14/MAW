[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/commit-activity)
[![GitHub contributors](https://img.shields.io/github/contributors/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/contributors/)
[![GitHub issues](https://img.shields.io/github/issues/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/issues/)

<p align="center"><img width="528" alt="MAW" src="https://user-images.githubusercontent.com/30716951/168855653-ae2efaa1-cbaf-4215-a04e-13bcd88ac46f.png"></p>


This repository hosts Metabolome Annotation Workflow (MAW). The workflow takes MS2 .mzML format data files as an input in R. It performs spectral database dereplication using R Package [Spectra](https://rformassspectrometry.github.io/Spectra/) and compound database dereplication using [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/). Final candidate selection is done in Python using [RDKit](https://www.rdkit.org/) and [PubChemPy](https://pubchempy.readthedocs.io/en/latest/).The classification of the tentative candidates from the input data are classified using [CANOPUS](https://bio.informatik.uni-jena.de/software/canopus/) and [ClassyFire](http://classyfire.wishartlab.com/).

## Install MAW

Install Docker on your [MAC OS](https://www.docker.com/get-started/) OR on [Linux](https://docs.docker.com/engine/install/ubuntu/). <br>

Since MAW is implemented in R and Python so we have two separate Docker images. To create a R-docker image, the following files are required which are present in MAW/R/cwl:
1. Dockerfile
2. Workflow_R_Functions.r
3. install_packages.R <br>

These files will create a docker image on your local system with the following command:
```shell
# clone the reporistory
git clone https://github.com/zmahnoor14/MAW.git
# go to the directory
cd MAW/R/cwl
# build the image
docker build -t mawr .
# go back to MAW
cd ..
cd ..
pwd # your local MAW Directory
/Users/my_name/Github/MAW
# enter the shell for the container
docker run --name container_name_of_your_choice -v $(pwd):/opt/workdir -i -t mawr /bin/bash
```


## Use Case

1. Download example .mzML file from this google drive link - [Link to MAW_ExampleData](https://drive.google.com/drive/folders/19huIu1sQ-bxgdx1C4nXZ75B_3eLpC9k-?usp=sharing). This file should be in the same folder as MAW/Tutorial.
2. Check the following files: MAW/Workflow_R_Workflow.r, MAW/Workflow_Python_Workflow.py
3. Check Tutorial files: MAW/Tutorial/Workflow_R_Script.r, MAW/Tutorial/Workflow_Python_Script.py
4. (optional to point 3)If you want to run Jupyter notebooks then check for following files: MAW/Tutorial/Workflow_R_Script.ipynb,  MAW/Tutorial/Workflow_Python_Script.ipynb
5. If following point 3, then MAW can simply be run using Rscript and python3 commands within Docker container. To do so, following commands are used:
```shell
# enter the shell for the container
docker run --name container_name_of_your_choice -v $(pwd):/opt/workdir -i -t mawr /bin/bash
#opt/workdir$
# for R
Rscript --no-save --no-restore --verbose Workflow_R_Script.r >outputFile.txt 2>&1 
# for Python
python3 Workflow_Python_Script.py
```
For details on using the workflow on Jupyter notebooks in a more interactive mode, please follow the Tutorial part on wiki page of this repository

## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
