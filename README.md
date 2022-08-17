[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/commit-activity)
[![GitHub contributors](https://img.shields.io/github/contributors/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/contributors/)
[![GitHub issues](https://img.shields.io/github/issues/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/issues/)

<p align="center"><img width="528" alt="MAW" src="https://user-images.githubusercontent.com/30716951/168855653-ae2efaa1-cbaf-4215-a04e-13bcd88ac46f.png"></p>

# Metabolome Annotation Workflow
 

This repository hosts Metabolome Annotation Workflow (MAW). The workflow has been developed using the LCMS-2 dataset from a marine diatom _Skeletonema marinoi_. The workflow takes .mzML format data files as an input in R and performs spectral database dereplication using R Package [Spectra](https://rformassspectrometry.github.io/Spectra/) and compound database dereplication using [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) and [MetFrag](https://ipb-halle.github.io/MetFrag/projects/metfragcl/) (with KEGG and PubChem). The results are saved as .csv files and are post processed in Python using [RDKit](https://www.rdkit.org/) and [PubChemPy](https://pubchempy.readthedocs.io/en/latest/).The classification of the tentative candidates from the input data are classified using [CANOPUS]() and [ClassyFire](http://classyfire.wishartlab.com/), with a python client [pybatchclassyfire](https://gitlab.unige.ch/Pierre-Marie.Allard/pybatchclassyfire/-/tree/master) for ClassyFire.

## Install MAW with Docker container
Install Docker on your MAC OS with (https://www.docker.com/get-started/) and for Linux with (https://docs.docker.com/engine/install/ubuntu/). <br>
Since MAW is implemented in R and Python so we have two separate Docker images.To create a R-docker image, the following files are required:
1. Dockerfile
2. Workflow_R_Functions.r
3. install_packages.R <br>
<<<<<<< HEAD
These files will create a docker image on your local system with the following command:
```shell
# build the image
docker build -t maw .
# run MAW in a jupyter notebook; change /workdir to your working directory
docker run -v /workdir:/workdir -i -t -p 8888:8888 maw /bin/bash -c "jupyter notebook --notebook-dir=/workdir --ip='*' --port=8888 --no-browser --allow-root"
=======
These files can be accessed via the folder R/cwl/.  A docker image will be created on your local system with these files using following command:
```shell
# clone the reporistory
# go to the directory
cd MAW/R/cwl
# build the image
docker build -t mawr .
# run MAW in a jupyter notebook; change /workdir to your working directory
docker run --name container_name_of_your_choice -v $(pwd):/opt/workdir -i -t mawr /bin/bash
>>>>>>> refs/remotes/origin/main
```
If you want to run the workflow in a docker container, then use the following command
```shell
docker run -i -t maw /bin/bash
#either run R or python3 within the container shell
```
If you want to use the working directory on your system rather than the docker container, then follow this command:
```shell
#create a directory named within the docker container
mkdir /mnt
#quit the session and type this command in your system shell
#check if you docker directory is your pwd and then run:
docker run -v $(pwd):/mnt -i -t maw /bin/bash
#either run R or python3 within the container shell
```
You are ready to use the workflow on a docker container on your system


## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
