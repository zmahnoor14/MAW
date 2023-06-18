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

#### Pull maw-r Docker Image
To pull MAW-R R-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-r:1.0.8 --platform=linux/amd64
```
This will creat a R-Docker image on your system. This image contains /opt/workdir as the working directory. This docker image will be used to run not only MAW-R but also MetFrag.

#### Pull maw-py Docker Image
To pull MAW-Py Python-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-py:1.0.7 --platform=linux/amd64
```
This will create a Python-Docker image on your system. This image contains /opt/workdir as the working directory.

#### Files

1. Download [cwl/Workflow_R_Script_all_MetFrag.r](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/Workflow_R_Script_all_MetFrag.r), and [cwl/Workflow_Python_Script_all_docker.py](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/Workflow_Python_Script_all_docker.py).
2. Download the databases gnps.rda, hmdb.rda, and mbankNIST.rda from https://zenodo.org/record/7519270. This submission contains GNPS saved at 2023-01-09 15:24:46, HMDB Current Version (5.0), MassBank version 2022.12.
3. Download COCONUT database to use with MetFrag from https://zenodo.org/record/7704937. If the user wants to use any other local libary, it is possible to use that instead of COCONUT. We recommend that the local file should be a csv file with atleast the following columns:     "Identifier"    "InChI"   "SMILES"    "molecular_weight". If you don't have information on all columns, these can be calculated with either RDKit or PubChempy automatically or can be done manually. If any important column is missing, MetFrag will give an error stating the name of the column.
4. Download MetFragCL [Jar file](https://github.com/ipb-halle/MetFragRelaunched/releases/download/v.2.5.0/MetFragCommandLine-2.5.0.jar). 
5. The most important file is the query input .mzML file

## Run MAW

### With CWL

It is important to keep all the CWL files in one folder. For other files you can provide a path wherever your files are.

```yaml
cwltool --cachedir cache maw.cwl /path/to/input.yaml --outputdir path/output/dir
```
maw.cwl will run the workflow taking the inputs from yaml file and store all the intermediate results in cache (for provenance and detailed results) and store the final results in outputdir. Repeat this for next .mzML file and so on.

### With Docker

#### Run MAW-R and MetFrag

```
docker run --name sample_name_maw-r -v $(pwd):/opt/workdir/ --platform linux/amd64  -it zmahnoor/maw-r:1.0.8 /bin/bash
```
Now you have entered a docker container with the name `sample_name_maw-r` and your pwd is mounted to /opt/workdir, so all the files necessary for docker will be mounted to /opt/workdir given that they were in your pwd. Note that using Workflow_R_Script_all_MetFrag.r you will run MAW-R and MetFrag in one script. pwd is your present working directory.
```
Rscript --no-save --no-restore --verbose Workflow_R_Script_all_MetFrag.r your_file_name.mzML gnps.rda hmdb.rda mbankNIST.rda give_an_id 15 TRUE coconut COCONUT_Jan2022.csv your_file_name/insilico/metparam_list.txt MetFragCommandLine-2.5.0.jar >script_output_Run_1.txt 2>&1 &

```
Only replace the `your_file_name` from the above command and it should work, also you can change the `give_an_id` and give your stdout file at the end a name of your choice. Also, if you use any local database, change the path to coconut and also change the `coconut` in the command line to any name you want to give to your database. The results will be stored like: /insilico/MetFrag/your_databasename.
<br>
one precursor mass takes 2 minutes on an Ubuntu system with 64GB RAM to run Workflow_R_Script_all_MetFrag.r. But since the files have more than one precursor mass, I recommend to run: `disown PID`. You will get a PID number when you enter and run the above Rscript command. The workflow will run on the background and you can exit your terminal, but before that you also need to exit the docker container as well. To do that, press CTRL+q and CTRL+p. Now the docker container runs in daemon mode.  To re enter the container:

```
docker exec -it sample_name_maw-r /bin/bash
```
The results will be saved in pwd.

#### Run MAW-Py
```
docker run --name sample_name_maw-py -v $(pwd):/opt/workdir/ --platform linux/amd64 -i -t zmahnoor/maw-py:1.0.7 /bin/bash
```
Once you enter the docker container, you can run the script and check the results.
Please make the following change in the script before running:

```
metfrag_candidate_list = pd.read_csv("ms2_spectra_ENDOpos/insilico/metparam_list.txt", sep = "\t", header=None, names=["metfrag_csv"])
```
TO
```
metfrag_candidate_list = pd.read_csv("your_file_name/insilico/metparam_list.txt", sep = "\t", header=None, names=["metfrag_csv"])
```
Then run the following command on terminal:
```
python3.10 Workflow_Python_Script_all_docker.py --file_id give_an_id --msp_file your_file_name/spectral_dereplication/spectral_results.csv --gnps_dir your_file_name/spectral_dereplication/GNPS --hmdb_dir your_file_name/spectral_dereplication/HMDB --mbank_dir your_file_name/spectral_dereplication/MassBank --ms1data your_file_name/insilico/MS1DATA.csv --score_thresh 0.75 > maw_py_log_file.txt 2>&1 &
```
This command will run the Workflow_Python_Script.py which will postprocess the results obtained earlier with MAW-R. In order to leave the container without killing the process, add & at the end of commands inside of the container. Then disown the PID and leave the container with CTRL+p and CTRL+q. <br>

The final results can be seen in the give_an_id_mergedResults-with-one-Candidates.csv.
<br>

<<<<<<< HEAD
For further information on the scripts used for the Bechmark datasets to validate results from MAW, please refer to the [MAW-Benchmark repository](https://github.com/zmahnoor14/MAW-Benchmark).

> **Note**
> For details on using the workflow on Jupyter notebooks in a more interactive mode, please follow the Tutorial part on wiki page of this repository.
=======

### Updated with SIRIUS5
It is possible to run SIRIUS and obtain results from [SIRIUS5](https://github.com/zmahnoor14/MAW/blob/provenance/cwl/Run_Sirius.r), however, MAW currently doesn't perform candidate selection with SIRIUS5 yet.
SIRIUS version 5 can be used using another docker container made for SIRIUS5. to run this, it is important that you already have .ms files in the SIRIUS folder generated from the MAW-R part of the workflow in the correct order of directory. Also, you would need to add your SIRIUS credentials after you enter the docker container. There is an R script called Run_SIRIUS5.r which contains the libraries, the function to run SIRIUS5 from R and a function call with your data. 

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

 SL means suspect list but this function is currently being updated by SIRIUS5 and shoulbe kept `FALSE`. profile can be the mass spectrometer used, either "orbitrap" or "qtof". db can be "ALL", "bio", "coconut" or any relevant database that is already provided by SIRIUS5. Also, for files, provide the full path of your './insilico/MS1DATA_SiriusP.tsv'. Run the function for each of your input files.

Enter docker container:
```
docker run --platform linux/x86_64 -i -t zmahnoor/maw-sirius5-old /bin/bash
```

Enter login info:

```
sirius login -u user@email.com -p --show
```
Th prompt will ask for your password, so enter the password and hit enter.

```
cd /data
Rscript Run_Sirius.r --files /your_file_name/insilico/MS1DATA_SiriusP.tsv --QC FALSE --SL FALSE SL_Path NA --profile orbitrap --db coconut >outputFile.txt 2>&1
```

> **Important Note**
> For details on using the workflow on Jupyter notebooks in a more interactive mode, please follow the Tutorial part on wiki page of this repository
> The current version only takes on .mzML file; the Batch processing of multiple .mzML files will be available soon
> At the moment, MAW runs MetFrag instead of SIRIUS5, however, MAW-R still generates input .ms files that can be used with SIRIUS5. We will very soon integrate SIRIUS5 in MAW (however, a SIRIUS account is required now which is free for academic use).
> Please add an issue for any bug report or any problems while running the workflow or write me an email at mahnoor.zulfiqar@uni-jena.de.
>>>>>>> provenance

## Citation
Zulfiqar, M., Gadelha, L., Steinbeck, C. et al. MAW: the reproducible Metabolome Annotation Workflow for untargeted tandem mass spectrometry. J Cheminform 15, 32 (2023). https://doi.org/10.1186/s13321-023-00695-y

## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
