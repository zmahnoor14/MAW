[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/commit-activity)
[![GitHub contributors](https://img.shields.io/github/contributors/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/contributors/)
[![GitHub issues](https://img.shields.io/github/issues/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/issues/)

<p align="center"><img width="528" alt="MAW" src="https://user-images.githubusercontent.com/30716951/168855653-ae2efaa1-cbaf-4215-a04e-13bcd88ac46f.png"></p>

# Metabolome Annotation Workflow
 

This repository hosts Metabolome Annotation Workflow (MAW). The workflow has been developed using the LCMS-2 dataset from a marine diatom _Skeletonema marinoi_. The workflow takes .mzML format data files as an input in R and performs spectral database dereplication using R Package [Spectra](https://rformassspectrometry.github.io/Spectra/) and compound database dereplication using [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) and [MetFrag](https://ipb-halle.github.io/MetFrag/projects/metfragcl/) (with KEGG and PubChem). The results are saved as .csv files and are post processed in Python using [RDKit](https://www.rdkit.org/) and [PubChemPy](https://pubchempy.readthedocs.io/en/latest/).The classification of the tentative candidates from the input data are classified using [CANOPUS]() and [ClassyFire](http://classyfire.wishartlab.com/), with a python client [pybatchclassyfire](https://gitlab.unige.ch/Pierre-Marie.Allard/pybatchclassyfire/-/tree/master) for ClassyFire.

## Install MAW with Docker container
Install Docker on your MAC OS with (https://www.docker.com/get-started/) and for Linux with (https://docs.docker.com/engine/install/ubuntu/). <br>
To create a docker image, the following files are required:
1. Dockerfile
2. maw.yml
3. install_packages.R <br>
These files will create a docker image on your local system with the following command:
```shell
# build the image
docker build -t maw .
# run MAW in a jupyter notebook; change /workdir to your working directory
docker run -v /workdir:/workdir -i -t -p 8888:8888 maw /bin/bash -c "jupyter notebook --notebook-dir=/workdir --ip='*' --port=8888 --no-browser --allow-root"
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

## Install MAW with Conda Environment
To create a conda environment, install conda from (https://www.anaconda.com/products/distribution). Once installed, type:
```shell
conda init # add to bashrc
conda
```
To install the packages for MAW, using the maw.yml environment file which is created using conda, follow the command:
```
conda env create -f maw.yml
activate mawRpy # which is the name of the environment in the maw.yml file
```
Currently, MsBackendMsp requires R 4.2, but to avoid the need to install R 4.2 which is unavailable on conda currently), use Docker image. The following code is only possible if you have an installation folder for MsBackednMsp saved (as it was previously available to download with 4.1.2). Follow the following instructions on the terminal to make the package available in your current R installation if you have a folder installed.
```shell
conda activate myenv
echo $CONDA_PREFIX
# you will receive a path, where you can keep the MsBackendMsp folder
```
There are some packages that are not available via any channel on conda OR there are errors downloading these packages with conda. To install such R packages, run the install_packages.R file within the mawRpy environment. For Python, install the following packages using ```pip3 install``` command.
```
pip3 install rdkit-pypi pubchempy requests_cache pybatchclassyfire
```

### Install SIRIUS
If the conda environment installation method is used, please install the latest SIRIUS version with <https://bio.informatik.uni-jena.de/software/sirius/> for Linux or MAC OS.
1. Installation with Linux
  ```shell
  echo $PATH # check which is already added to PATH variable
  sudo emacs /etc/environment
  PATH="usr/s_cost/sirius-gui/bin"
  source /etc/environment
  ```
2. Installation with MAC
  To use this command line, add the following path to either .bash_profile or .zprofile. You can find these files easily using FileZilla <https://filezilla-project.org/>. open your user name folder. Open your .bash_profile, add this to the $PATH variable.
  ```shell
  PATH="/usr/s_cost/sirius.app/Contents/MacOS/:${PATH}"
  export PATH
  ```
The path may differ, so check your Sirius installation folder to get the correct path name.

## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
