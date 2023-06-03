[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/commit-activity)
[![GitHub contributors](https://img.shields.io/github/contributors/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/contributors/)
[![GitHub issues](https://img.shields.io/github/issues/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/issues/)
[![DOI](https://zenodo.org/badge/438345970.svg)](https://zenodo.org/badge/latestdoi/438345970)

<p align="center"><img width="528" alt="MAW" src="https://user-images.githubusercontent.com/30716951/168855653-ae2efaa1-cbaf-4215-a04e-13bcd88ac46f.png"></p>


This repository hosts Metabolome Annotation Workflow (MAW). The workflow takes MS2 .mzML format data files as an input in R. It performs spectral database dereplication using R Package [Spectra](https://rformassspectrometry.github.io/Spectra/) and compound database dereplication using [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) and [MetFrag](https://ipb-halle.github.io/MetFrag/) . Final candidate selection is done in Python using [RDKit](https://www.rdkit.org/) and [PubChemPy](https://pubchempy.readthedocs.io/en/latest/).The classification of the tentative candidates from the input data are classified using [CANOPUS](https://bio.informatik.uni-jena.de/software/canopus/) and [ClassyFire](http://classyfire.wishartlab.com/).


## Installation

### Common Workflow Language (CWL)

Users can run MAW with CWL with just one command without pulling individual docker images. 
1. install [cwltool](https://github.com/common-workflow-language/cwltool#install)
2. Download [cwl/Workflow_R_Script_all.r](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/Workflow_R_Script_all.r) and [cwl/Workflow_Python_Script_all.py](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/Workflow_Python_Script_all.py). These will be the scripts used by MAW CWL.
3. Download the databases gnps.rda, hmdb.rda, and mbankNIST.rda from https://zenodo.org/record/7519270. This submission contains GNPS saved at 2023-01-09 15:24:46, HMDB Current Version (5.0), MassBank version 2022.12.
4. Download COCONUT database to use with MetFrag from https://zenodo.org/record/7704937. If the user wants to use any other local libary, it is possible to use that instead of COCONUT. We recommend that the local file should be a csv file with atleast the following columns:     "Identifier"    "InChI"   "SMILES"    "molecular_weight". If you don't have information on all columns, these can be calculated with either RDKit or PubChempy automatically or can be done manually. If any important column is missing, MetFrag will give an error stating the name of the column.
5. The most important file is the query input .mzML file
6. Download [cwl/maw-r.cwl](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/maw-r.cwl), [cwl/maw-py.cwl](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/maw-py.cwl), [cwl/maw-metfrag.cwl](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/maw-metfrag.cwl), [cwl/maw.cwl](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/maw.cwl)
7. To make an input yaml, take an example from the [cwl/template_maw.yml](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/template_maw.yaml)


### Docker containers

Install Docker on your [MAC OS](https://www.docker.com/get-started/) OR on [Linux](https://docs.docker.com/engine/install/ubuntu/). <br>

At the moment MAW runs R script, MetFrag, and Python Script. Note that SIRIUS5 is not yet implemented yet, but is in progress.
So, here we have two docker containers and 3 scripts --> (MAW-R, MetFrag.r, and MAW-Py).

### Pull maw-r Docker Image
To pull MAW-R R-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-r:1.0.8 --platform=linux/amd64
```
This will creat a R-Docker image on your system. This image contains /opt/workdir as the working directory. This docker image will be used to run not only MAW-R but also MetFrag.

### Pull maw-py Docker Image
To pull MAW-Py Python-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-py:1.0.7 --platform=linux/amd64
```
This will create a Python-Docker image on your system. This image contains /opt/workdir as the working directory.

## Run MAW

### With CWL

It is important to keep all the CWL files in one folder. For other files you can provide a path wherever your files are.

```yaml
cwltool --cachedir cache maw.cwl /path/to/input.yaml --outputdir path/output/dir
```
maw.cwl will run the workflow taking the inputs from yaml file and store all the intermediate results in cache (for provenance and detailed results) and store the final results in outputdir. Repeat this for next .mzML file and so on.

### With Docker

Docker is a bit more complicated to use but gives more freedom and options right now as cwl is still in its development phase.

```
docker run --name example_maw-r -i -t zmahnoor/maw-r:1.0.0 /bin/bash
```
Once you enter the docker container, you can run the script and check the results.
```
Rscript --no-save --no-restore --verbose Workflow_R_Script.r >outputFile.txt 2>&1
```
This command will run the Workflow_R_Script.r which is an example script for /data/example_Tyrosine.mzML. The calculation takes about 2 minutes on an Ubuntu system with 64GB RAM.

### [NEW] 2. Updated with SIRIUS5

SIRIUS version 5 can be used using another docker container made for SIRIUS5. to run this, it is important athat you alsready have .md files in the SIRIUS folder generated from the MAW-R part of the workflow in the correct order of directory. Also, you would need to add your SIRIUS credentials after you enter the docker container. THere is an R script called Run_SIRIUS5.r which contains the libraries, the function to run SIRIUS5 from R and a function call with your data. 

```R
run_sirius(files= './insilico/MS1DATA_SiriusP.tsv',
           ppm_max = 5, 
           ppm_max_ms2 = 15, 
           QC = FALSE, 
           SL = FALSE, 
           SL_path = NA, 
           candidates = 30, 
           profile = "orbitrap", 
           db = "coconut")
```

Please check the R script and make the changes in the parameters above according to your data needs. SL means suspect list but this function is currently being updated by SIRIUS5 and shoulbe kept ```FALSE```. profile can be the mass spectrometer used, either "orbitrap" or "qtof". db can be "ALL", "bio", "coconut" or any relevant database that is already provided by SIRIUS5. Also, for files, provide the full path of your './insilico/MS1DATA_SiriusP.tsv'. Run the function for each of your input files.

Enter docker container:
```
docker run --platform linux/x86_64 -i -t zmahnoor/maw-sirius5-old /bin/bash
```

Enter login info:

```
sirius login -u user@email.com -p --show
```
Th prompt will ask for your password, so enter the pssowrd and hit enter.

```
cd /data
Rscript --no-save --no-restore --verbose Run_SIRIUS5.r >outputFile.txt 2>&1
```

### 3. Run maw-py
```
docker run --name example_maw-py -i -t zmahnoor/maw-py:1.0.0 /bin/bash
```
Once you enter the docker container, you can run the script and check the results.
```
python3 Workflow_Python_Script.py
```
This command will run the Workflow_Python_Script.py which is an example script for /data/Example_Tyrosine/ obtained earlier with MAW-R. In order to leave the container without killing the process, add $ at the end of commands inside of the container. Then disown the PID and leave the container with CTRL+P and CTRL+Q. <br>

The results can be seen in the /data/final_candidates.csv.
<br>

> **Important Note**
> For details on using the workflow on Jupyter notebooks in a more interactive mode, please follow the Tutorial part on wiki page of this repository
> The current version only takes on .mzML file; the Batch processing of multiple .mzML files will be available soon
> At the moment, MAW runs MetFrag instead of SIRIUS5, however, MAW-R still generates input .ms files that can be used with SIRIUS5. We will very soon integrate SIRIUS5 in MAW (however, a SIRIUS account is required now which is free for academic use).

## Citation
Zulfiqar, M., Gadelha, L., Steinbeck, C., Sorokina, M., & Peters, K. (2022). Metabolome Annotation Workflow (MAW) (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.7148450

## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
