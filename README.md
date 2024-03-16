[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/commit-activity)
[![GitHub contributors](https://img.shields.io/github/contributors/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/contributors/)
[![GitHub issues](https://img.shields.io/github/issues/zmahnoor14/MAW.svg)](https://GitHub.com/zmahnoor14/MAW/issues/)
[![DOI](https://zenodo.org/badge/438345970.svg)](https://zenodo.org/badge/latestdoi/438345970)
[![FAIR checklist badge](https://fairsoftwarechecklist.net/badge.svg)](https://fairsoftwarechecklist.net/v0.2?f=31&a=32113&i=31011&r=133)

<p align="center"><img width="528" alt="MAW" src="https://user-images.githubusercontent.com/30716951/168855653-ae2efaa1-cbaf-4215-a04e-13bcd88ac46f.png"></p>

This repository hosts Metabolome Annotation Workflow (MAW). The workflow takes MS2 .mzML format data files as an input in R. It performs spectral database dereplication using R Package [Spectra](https://rformassspectrometry.github.io/Spectra/) and compound database dereplication using [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) OR [MetFrag](https://ipb-halle.github.io/MetFrag/). Final candidate selection is done in Python using [RDKit](https://www.rdkit.org/) and [PubChemPy](https://pubchempy.readthedocs.io/en/latest/). The classification of the tentative candidates from the input data are classified using [CANOPUS](https://bio.informatik.uni-jena.de/software/canopus/) and [ClassyFire](http://classyfire.wishartlab.com/), which will be updated to accomodate [NPClassifier](https://npclassifier.ucsd.edu/) soon.


## Installation

MAW is available as docker containers and as CWL workflow description. It can be used via both channels, however at the moment, SIRIUS is only accommodated with the docker containers, (the workflow will be completely operable in CWL in the future). At the moment, the CWL version can be used with MetFrag. Given below are the instructions on installation either via CWL (option to use MetFrag with all the additional features of the workflow) or via Docker (options to run either SIRIUS or MetFrag in addition to all other features of the workflow). With CWL, multiple files can be executed with MAW in parallel. it is recommended to use HPC for bigger files via a batch system such as SLURM (>= 10 precursor masses).

### Common Workflow Language (CWL)

Users can run MAW with CWL with just one command without pulling individual docker images. 
1. install [cwltool](https://github.com/common-workflow-language/cwltool#install)
2. Download [cwl/Workflow_R_Script_all.r](https://github.com/zmahnoor14/MAW/blob/main/cwl/Workflow_R_Script_all.r) and [cwl/Workflow_Python_Script_all.py](https://github.com/zmahnoor14/MAW/blob/main/cwl/Workflow_Python_Script_all.py). These will be the scripts used by MAW CWL.
3. Download the databases gnps.rda, hmdb.rda, and mbankNIST.rda from https://zenodo.org/record/7519270. This submission contains GNPS saved at 2023-01-09 15:24:46, HMDB Current Version (5.0), MassBank version 2022.12.
4. Download COCONUT database to use with MetFrag from https://zenodo.org/record/7704937. If the user wants to use any other local library, it is possible to use that instead of COCONUT. We recommend that the local file should be a csv file with atleast the following columns:     "Identifier"    "InChI"   "SMILES"    "molecular_weight". If you don't have information on all columns, these can be calculated with either RDKit or PubChempy automatically or can be done manually. If any important column is missing, MetFrag will give an error stating the name of the missing column.
5. The most important file is the query input .mzML files
6. Download [cwl/maw-r.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-r.cwl), [cwl/maw-py.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-py.cwl), [cwl/maw-metfrag.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-metfrag.cwl), [cwl/maw_single.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw_single.cwl), and [cwl/maw.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw.cwl)
7. To make an input yaml, take an example from the [cwl/maw-inputs-multiple_ara.yaml](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-inputs-multiple_ara.yaml)


### Docker containers

Install Docker on your [MAC OS](https://www.docker.com/get-started/) OR on [Linux](https://docs.docker.com/engine/install/ubuntu/). <br>

#### MAW-with-MetFrag

* Pull maw-r Docker Image
To pull MAW-R R-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-r:1.0.8 --platform=linux/amd64
```
This will create a R-Docker image on your system. This image contains /opt/workdir as the working directory. This docker image will be used to run not only MAW-R but also MetFrag.

* Pull maw-py Docker Image
To pull MAW-Py Python-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-py:1.0.7 --platform=linux/amd64
```
This will create a Python-Docker image on your system. This image contains /opt/workdir as the working directory.

* Files

1. Download [Docker/MetFrag/Workflow_R_Script_all_MetFrag.r](https://github.com/zmahnoor14/MAW/blob/main/Docker/MetFrag/Workflow_R_Script_all_MetFrag.r), and [Docker/MetFrag/Workflow_Python_Script_all_MetFrag.py](https://github.com/zmahnoor14/MAW/blob/main/Docker/MetFrag/Workflow_Python_Script_all_docker.py).
2. Download the databases gnps.rda, hmdb.rda, and mbankNIST.rda from https://zenodo.org/record/7519270. This submission contains GNPS saved at 2023-01-09 15:24:46, HMDB Current Version (5.0), MassBank version 2022.12.
3. Download COCONUT database to use with MetFrag from https://zenodo.org/record/7704937. If the user wants to use any other local library, it is possible to use that instead of COCONUT. We recommend that the local file should be a csv file with atleast the following columns:     "Identifier"    "InChI"   "SMILES"    "molecular_weight". If you don't have information on all columns, these can be calculated with either RDKit or PubChempy automatically or can be done manually. If any important column is missing, MetFrag will give an error stating the name of the column.
4. Download MetFragCL [Jar file](https://github.com/ipb-halle/MetFragRelaunched/releases/download/v.2.5.0/MetFragCommandLine-2.5.0.jar). 
5. The most important file is the query input .mzML file


#### MAW-with-SIRIUS
* Pull maw-r Docker Image
To pull MAW-R R-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-r:1.0.8 --platform=linux/amd64
```
* Pull maw-sirius Docker Image
```
docker pull zmahnoor/maw-sirius5-old:1.0.0 --platform=linux/amd64
```

* Pull maw-py Docker Image
To pull MAW-Py Python-docker image, run the following command on your terminal:
```
docker pull zmahnoor/maw-py:1.0.7 --platform=linux/amd64
```

* Files

1. Download [Docker/SIRIUS5/Workflow_R_Script_all_SIRIUS.r](https://github.com/zmahnoor14/MAW/blob/main/Docker/SIRIUS5/Workflow_R_Script_all_SIRIUS.r), and [Docker/SIRIUS5/Workflow_Python_Script_all_SIRIUS.py](https://github.com/zmahnoor14/MAW/blob/main/Docker/SIRIUS5/Workflow_Python_Script_all_SIRIUS.py).
2. Download the databases gnps.rda, hmdb.rda, and mbankNIST.rda from https://zenodo.org/record/7519270. This submission contains GNPS saved at 2023-01-09 15:24:46, HMDB Current Version (5.0), MassBank version 2022.12.
3. The most important file is the query input .mzML file

## Run MAW

### With CWL

It is important to keep all the CWL files (maw-r.cwl, maw-py.cwl, maw-metfrag.cwl, maw_single.cwl, and maw.cwl) in one folder. For other files, you can provide a path wherever your files are. The input yaml file is inspired by the [cwl/maw-inputs-multiple_ara.yaml](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-inputs-multiple_ara.yaml). Remember with CWL, only MetFrag version of MAW is available at the moment.

```yaml
cwltool --cachedir cache maw.cwl /path/input/maw-inputs-multiple_ara.yaml --outputdir path/output/dir
```
maw.cwl will run the workflow taking the inputs from yaml file store all the intermediate results in the cache (for provenance and detailed results) and store the final results in outputdir. This command will run all the input mzML files in a loop. To run in parallel, ```--parallel``` is added to the above command after ```cwltool```.

### With Docker

#### Run MAW (with MetFrag option)
* Run MAW-R. The MAW-R uses Workflow_R_Script_all_MetFrag.r, which not only performs spectral database dereplication but also runs MetFrag.
```
docker run --name sample_name_maw-r -v $(pwd):/opt/workdir/data --platform linux/amd64  -it zmahnoor/maw-r:1.0.8 /bin/bash
```
Now you have entered a docker container with the name `sample_name_maw-r` and your pwd is mounted to /opt/workdir/data, so all the files necessary for docker will be mounted to /opt/workdir/data given that they were in your pwd. pwd is your present working directory. Follow the code given below once you enter the docker container from the above command.
```
cd ..
cd data
Rscript --no-save --no-restore --verbose Workflow_R_Script_all_MetFrag.r your_file_name.mzML gnps.rda hmdb.rda mbankNIST.rda 15 TRUE coconut COCONUT_Jan2022.csv your_file_name/insilico/metparam_list.txt MetFragCommandLine-2.5.0.jar >script_output_Run_1.txt 2>&1 &

```
Workflow_R_Script_all_MetFrag.r = script<br>
your_file_name.mzML = input file <br>
gnps.rda = gnps database R object <br>
hmdb.rda = hmdb database R object <br>
mbankNIST.rda = massbank database R object <br>
15 = ppm for spectral database dereplication precursor matches <br>
TRUE = collision energy <br>
coconut = id for database used <br>
COCONUT_Jan2022.csv = path to coconut database as CSV <br>
your_file_name/insilico/metparam_list.txt = an output file generated during MAW-R runs. This file will be used to run MetFrag at the end. You will only get this file after the second function runs. <br>
MetFragCommandLine-2.5.0.jar = path to the metfrag jar file, used for running metfrag.<br>

Only replace the `your_file_name` from the above command and it should work and give your stdout file at the end a name of your choice. Also, if you use any local database, change the path to coconut and also change the `coconut` in the command line to any name you want to give to your database. The results will be stored like: /insilico/MetFrag/your_databasename.
<br>
one precursor mass takes 2 minutes on an Ubuntu system with 64GB RAM to run Workflow_R_Script_all_MetFrag.r. But since the files have more than one precursor mass, I recommend running: `disown PID`. You will get a PID number when you enter and run the above Rscript command. The workflow will run in the background and you can exit your terminal, but before that, you also need to exit the docker container as well. To do that, press CTRL+q and CTRL+p. Now the docker container runs in daemon mode.  To re-enter the container:

```
docker exec -it sample_name_maw-r /bin/bash
```
The results will be saved in pwd.

* Run MAW-Py. First run the docker container for MAW-Py
```
docker run --name sample_name_maw-py -v $(pwd):/opt/workdir/data --platform linux/amd64 -i -t zmahnoor/maw-py:1.0.7 /bin/bash
```
Once you enter the docker container, you can run the script and check the results. Run the following command

```
cd data
python3.10 Workflow_Python_Script_all_MetFrag.py --msp_file your_file_name/spectral_dereplication/spectral_results.csv --gnps_dir your_file_name/spectral_dereplication/GNPS --hmdb_dir your_file_name/spectral_dereplication/HMDB --mbank_dir your_file_name/spectral_dereplication/MassBank --ms1data your_file_name/insilico/MS1DATA.csv --score_thresh 0.75 > maw_py_log_file.txt 2>&1 &
```
'--msp_file', type=str, help='path to spec result CSV file' <br>
'--gnps_dir', type=str, help='path to GNPS directory' <br>
'--hmdb_dir', type=str, help='path to HMDB directory' <br>
'--mbank_dir', type=str, help='path to MassBank directory' <br>
'--ms1data', type=str, help='path to MS1 data CSV file' <br>
'--score_thresh', type=float, default=0.75, help='score threshold for MetFrag results (default: 0.75)' <br>
msp_file, gnps_dir,hmdb_dir,mbank_dir, ms1data are generated by MAW-R. The above command gives the paths to the results generated. e.g.: only replace your_file_name with your mzml file name, without the .mzML extension. <br>
This command will run the Workflow_Python_Script.py which will postprocess the results obtained earlier with MAW-R. To leave the container without killing the process, add & at the end of commands inside of the container. Then disown the PID and leave the container with CTRL+p and CTRL+q. <br>

The final results can be seen in the file_id_mergedResults-with-one-Candidates.csv.
<br>

### Run MAW (with SIRIUS5 option)
So with SIRIUS, there are three docker containers and three steps: 1) Run MAW-R script from MAW/Docker/SIRIUS5 within the docker container, 2) Run the docker container for SIRIUS, with the instructions given in this section, 3) Run MAW-Py script from MAW/Docker/SIRIUS5 within the docker container. The instructions to run SIRIUS, are given below.

* Run MAW_R, using the file [Docker/SIRIUS5/Workflow_R_Script_all_SIRIUS.r](https://github.com/zmahnoor14/MAW/blob/main/Docker/SIRIUS5/Workflow_R_Script_all_SIRIUS.r):
```
docker run --name sample_name_maw-r -v $(pwd):/opt/workdir/data --platform linux/amd64  -it zmahnoor/maw-r:1.0.8 /bin/bash
```
once you enter the docker container, add the following commands: 
```
cd ..
cd data
Rscript --no-save --no-restore --verbose Workflow_R_Script_all_MetFrag.r your_file_name.mzML gnps.rda hmdb.rda mbankNIST.rda 15 TRUE >script_output_Run_1.txt 2>&1 &

```
Workflow_R_Script_all_MetFrag.r = script<br>
your_file_name.mzML = input file <br>
gnps.rda = gnps database R object <br>
hmdb.rda = hmdb database R object <br>
mbankNIST.rda = massbank database R object <br>
15 = ppm for spectral database dereplication precursor matches <br>
TRUE = collision energy <br>

* Run SIRIUS5 <br>
It is possible to run SIRIUS and obtain results from [SIRIUS5](https://github.com/zmahnoor14/MAW/blob/main/Docker/SIRIUS5/Run_Sirius.r), SIRIUS version 5 can be used using another docker container made for SIRIUS5. To run this, it is important that you already have .ms files in the SIRIUS folder generated from the MAW-R part of the workflow in the correct order of directory --> "file.mzML/insilico/SIRIUS/no_isotope (because so far no isotopic information is included). Also, you would need to add your SIRIUS credentials after you enter the docker container. There is an R script called Run_SIRIUS5.r which contains the libraries, the function to run SIRIUS5 from R, and a function call with your data. This file is already in the docker container. So following the steps, it is possible to run SIRIUS5.

1. Run a container with the working directory where you file.mzML and its result directory /file are stored:
```
docker run --name name_your_container -v $(pwd):/opt/workdir/data --platform linux/amd64 -it zmahnoor/maw-sirius5-old:1.0.0 /bin/bash
```
2. Put your credentials, this is important whenever you create a new container.
Enter login info:

```
sirius login -u user@email.com -p --show
```
The prompt will ask for your password, so enter the password and hit enter.

3. When you enter the container, the current directory is /opt/workdir/data. An example of the SIRIUS input files is: /opt/workdir/data/file/insilico/SIRIUS/no_isotope/example.ms<br>
```
cd data
Rscript Run_Sirius.r /opt/workdir/data/your_file_name/insilico/MS1DATA_SiriusP.tsv orbitrap coconut >outputFile.txt 2>&1&
```
The arguments after the R script can be explained by the R function that runs SIRIUS within R.
```R
run_sirius(files= '/opt/workdir/data/your_file_name/insilico/MS1DATA_SiriusP.tsv',
           ppm_max = 5, 
           ppm_max_ms2 = 15, 
           QC = FALSE, 
           SL = FALSE, 
           SL_path = NA, 
           candidates = 30, 
           profile = "orbitrap", 
           db = "coconut")
```

The first user-defined argument is named files, which is generally named like this in the workflow: /opt/workdir/data/your_file_name/insilico/MS1DATA_SiriusP.tsv and is generated during the MAW_R script. It contains the .ms input files and their paths, and the corresponding .json output directory where SIRIUS writes its outputs. The next argument is QC which remains `FALSE`, the SL means suspect list but this function is currently being updated by SIRIUS5 and should be kept `FALSE`, and hence the SL_path is NA. Profile can be the mass spectrometer used, either "orbitrap" or "qtof". db can be "ALL", "bio", "coconut" or any relevant database that is already provided by SIRIUS5. 

* Run MAW-Py<br>
```
docker run --name sample_name_maw-py -v $(pwd):/opt/workdir/data --platform linux/amd64 -i -t zmahnoor/maw-py:1.0.7 /bin/bash
```
Once you enter the docker container, you can run the script and check the results. Run the following command

```
cd data
python3.10 Workflow_Python_Script_all_SIRIUS.py --msp_file your_file_name/spectral_dereplication/spectral_results.csv --gnps_dir your_file_name/spectral_dereplication/GNPS --hmdb_dir your_file_name/spectral_dereplication/HMDB --mbank_dir your_file_name/spectral_dereplication/MassBank --ms1data your_file_name/insilico/MS1DATA.csv --db coconut > maw_py_log_file.txt 2>&1 &
```
'--msp_file', type=str, help='path to spec result CSV file' <br>
'--gnps_dir', type=str, help='path to GNPS directory' <br>
'--hmdb_dir', type=str, help='path to HMDB directory' <br>
'--mbank_dir', type=str, help='path to MassBank directory' <br>
'--ms1data', type=str, help='path to MS1 data CSV file' <br>
'--db', type=str, default='coconut', help='whichever database was used in SIRIUS, mostly coconut database' <br>

> **Important Note**
> The current version only takes on .mzML files; the Batch processing of multiple .mzML files is only available with CWL --parallel mode, which is only recommended with a batching system such as slurm on HPC.
> Please add an issue for any bug report or any problems while running the workflow or write me an email at mahnoor.zulfiqar@uni-jena.de OR zmahnoor14@gmail.com.

## Citation
Zulfiqar, M., Gadelha, L., Steinbeck, C. et al. MAW: the reproducible Metabolome Annotation Workflow for untargeted tandem mass spectrometry. J Cheminform 15, 32 (2023). https://doi.org/10.1186/s13321-023-00695-y

## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
