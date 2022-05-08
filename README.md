[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zmahnoor14/MAW/graphs/commit-activity)

# Metabolome Annotation Workflow
 - (tested on Linux and MAC OS)

This repository hosts Metabolome Annotation Workflow (MAW). The workflow has been developed using the LCMS-2 dataset from a marine diatom _Skeletonema marinoi_. The workflow takes .mzML format data files as an input in R and performs spectral database dereplication using R Package [Spectra](https://rformassspectrometry.github.io/Spectra/) and compound database dereplication using [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) and [MetFrag](https://ipb-halle.github.io/MetFrag/projects/metfragcl/) (with KEGG and PubChem). The results are savd as .csv files and are post processed in Python using [RDKit](https://www.rdkit.org/) and [PubChemPy](https://pubchempy.readthedocs.io/en/latest/).The classification of the tentative candidates from the input data are classified using [CANOPUS]() and [ClassyFire](http://classyfire.wishartlab.com/), with a python client [pybatchclassyfire](https://gitlab.unige.ch/Pierre-Marie.Allard/pybatchclassyfire/-/tree/master) for ClassyFire.

## Install MAW with Docker container
Install Docker on your MAC OS with (https://www.docker.com/get-started/) and for Linux with (https://docs.docker.com/engine/install/ubuntu/). <br>
To create a docker image, following files are required:
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
Currently, MsBackendMsp requires R 4.2, but to avoid the need to install R 4.2 which is unavailable on conda currently), use Docker image. The following code is only possible if you have an installation folder for MsBackednMsp saved (as it was previously available to doenaload with 4.1.2). Follow the following instructions on the terminal to make the package available in your current R installation if you have a folder installed.
```shell
conda activate myenv
echo $CONDA_PREFIX
# you will receive a path, where you can keep the MsBackendMsp folder
```
There are some packages that are not available via any channel on conda OR there are errors downloading these packags with conda. To install such R packages, run the install_packages.R file within the mawRpy environment. For Python, install following packages using ```pip3 install``` command.
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
The path may be different, so check your Sirius installation folder to get the correct path name.

## Input files and Directories

An input directory (/input_dir) should have the following files.
1. All LCMS-2 spectra .mzML files
2. SIRIUS installation zip folder (only if you install via conda environment; no need to have this with docker image installation)
3. MetFrag jar file which can be downloaded from <https://github.com/ipb-halle/MetFragRelaunched/releases>.
4. MetFrag_AdductTypes.csv can be downloaded from <https://github.com/schymane/ReSOLUTION/blob/master/inst/extdata/MetFrag_AdductTypes.csv>
5. (optional) /QC folder containing all (MS1) QC .mzML files 
6. (optional) Suspect List in .csv format (important column - "SMILES")


## Functions to store online Spectral Databases to your R environment [R]
In order to download spectral databases, use the following function:
```
download_specDB(input_dir, db = "all")
```
This function will take a lot of computational resources. However, to skip this function, you can download the current versions of these databases from <https://zenodo.org/deposit/6528931> (linked embargoed). The databases are stored in the same input directory and in .rda format, as an R object.

## Functions to use a Suspect List(in-house library) [Python]
A suspect list of compounds can be used within the workflow to provide confidence to the predictions. This library is matched against results from Spectral and Compound Databases. MAW provides two functions to generate the input suspect list compounds for SIRIUS and MetFrag. The only important information in the SuspectList.csv file should be a column with SMILES. 
1. SIRIUS requires a folder with many .tpt files which contain fragmentation tree for each SMILES. To generate suspect list input for SIRIUS, use the function ```slist_sirius```. This function generates a result folder /input_dir/SL_Frag. 
```
slist_sirius(input_dir, slist_csv, substring = ["NA+", "Zn+"])
```
2. MetFrag requires a txt file with InChIKeys. This can be obtained with the function ```slist_metfrag```. 
```
slist_metfrag(input_dir, slist_csv, name)
```
## Tutorial of Workflow

### R
Follow the jupyter notebook: Workflow_R_Script.ipynb
1. Load Dependencies:
```R
library(Spectra)
library(MsBackendMgf)
library(MsBackendHmdb)
library(MsCoreUtils)
library(MsBackendMsp)
library(readr)
library(dplyr)
library(rvest)
library(stringr)
library(xml2)
```

2. Define the input directory. Make sure that you have input LCMS-2 spectra files in .mzML format in the input directory.
```R
input_dir <- "/usr/input_dir/"
```

3. Load the database .rda objects to the current R session.
```R
load(file = paste(input_dir,"gnps.rda", sep = ""))
load(file = paste(input_dir,"hmdb.rda", sep = ""))
load(file = paste(input_dir,"mbank_NIST.rda", sep = ""))
```
4. Load the functions file
```R
source(file = paste(input_dir, "Workflow_R_Functions.r", sep = ''))
```
5. Create a table that lists all the input .mzML files and the result directories (created with the same function) with the same name as the input file. The function ```ms2_rfilename``` gives an id to the file as well. 
```R
input_table <- data.frame(ms2_rfilename(input_dir))
```

6. (optional) If you are using QC files, follow these steps as well. QC samples are used to normalize the signals across all samples. So, generally, QC files contain all the signals from MS1 files with higher m/z resolution. These files also contain the isotopic peaks and can be used for formula identification in SIRIUS. 
    1. If you have files with positive and negative modes in one file, follow the first section of the code. It takes all files in the QC folder with a certain pattern (choose a pattern from your file names, could be even .mzML as this is the format of all files in the QC folder) and divides the pos and neg modes from each file and generates different mode .mzML files (01pos.mzmL and 01neg.mzML). These files are read by [CAMERA](https://www.bioconductor.org/packages/release/bioc/html/CAMERA.html), which is loaded within the function. The outputs are several .csv files with CAMERA results, which are merged into one .csv file for each model.
    2. If you have files with positive and negative modes in separate files, follow the second section of code. It takes one file in the QC folder with either pos or neg mode. These files are read by "CAMERA", which is loaded within the function.          
```R
# first section
cam_funcMode(path = paste(input_dir, "QC", sep =""), pattern = ".mzML")
merge_qc(path = paste(input_dir, "QC", sep =""))
```

```R
# second section
cam_func(path = "QC/", f = "01pos.mzML", mode = "pos")
cam_func(path = "QC/", f = "02neg.mzML", mode = "neg")
```

7. (optional, if followed 6.) Add the QC files to input_table. This can be done simply by adding one more column to the input_table and adding the pos files against pos mode LCMS-2 data and neg files against neg mode LCMS-2 files. An example is given below:             
```R
for (i in 1:nrow(input_table)){
    if (grepl("pos", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_pos.csv"
    }
    if (grepl("neg", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_neg.csv"
    }
}
```
  
8. After the previous steps, we have all the inputs, their directories and optionally the QC .csv files. Next is to initiate the workflow, (assuming we want to process only one LCMS-2 .mzML file). Use the spec_Processing function to read and pre-process the MS2 spectra. The output is processed spectra and a list of precursor m/z(s) present in the .mzML file.                                     
```R
spec_pr <- spec_Processing(as.character(input_table[i, "mzml_files"]), input_table[i, "ResultFileNames"])
```
  
9. Using the above-processed spectra, the workflow is ready to perform spectral database deprelication. In this workflow, "GNPS", "HMDB" and "MassBank" or all of these databases can be used. This function can be performed for all precursor m/z(s) in a loop. The function ```spec_dereplication``` takes one precursor m/z at a time, the result directory path, either from the input_table or as a whole written path e.g: "usr/s_cost/file_pos_01". The file_id is also taken from the input_table, but is customizable as any other id given as a string e.g: "IDfile_pos_01". ppmx is the ppm value used by the function to match the two spectra which have their m/z values at most 15 ppm apart to be considered a match. This results in a directory called spectral_dereplication and contains results in .csv files for each database.
```R
df_derep <- spec_dereplication(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "./"), sep =""), 
                                   proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "./"), sep =""),
                                   db = "all", 
                                   result_dir = input_table[i, "ResultFileNames"],
                                   file_id = input_table[i, "File_id"], 
                                   input_dir, 
                                   ppmx = 15)
```
                
10. In order to perform dereplication using compound databases, the Workflow uses SIRIUS and MetFrag. For these tools, we need MS-2 fragmentation peak lists. The next function ```ms2_peaks``` extracts MS-2 fragmentation spectra and stores the peaks for each precursor m/z in "/insilico/peakfiles_ms2/Peaks_01.txt".             
```R
spec_pr2 <- ms2_peaks(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "./"), sep =""), 
                          proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "./"), sep =""),
                          input_dir,
                          input_table[i, "ResultFileNames"],
                         file_id = input_table[i, "File_id"])
```
  
11. For Formula identification in SIRIUS, a fragmentation tree generates a tree score and isotopic peaks generate an isotopic score, both of which form the Sirius score. In order to utilize the isotopic score, isotopic peak annotation is required. For this purpose, the CAMERA results, are utilized here to extract the isotopic peaks for each precursor m/z and store them in "/insilico/peakfiles_ms2/Peaks_01.txt". The function ```ms1_peaks``` takes results from ms2_peaks, the associated QC CSV file, the result directory and whether the QC file is being used or not. If there is not QC .csv file, enter NA and QC = FALSE.
  
```R
ms1p <- ms1_peaks(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS2DATA.csv', sep = ""), "./"), sep =""), 
                      y = input_table[i, "qcCAM_csv"], 
                      input_table[i, "ResultFileNames"], 
                      input_dir, 
                      QC = FALSE)
```
        
        
11. To create SIRIUS input files, use ```sirius_param_files``` function, which takes the output table from ```ms1_peaks``` and the result directory. If a suspect list is used, then SL = TRUE. 
```R
sirius_param_files <- sirius_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA.csv', sep = ""), "./"), sep =""), 
                                       result_dir = input_table[i, 'ResultFileNames'], 
                                       input_dir,
                                       SL = TRUE)
```
  
12. Run SIRIUS. Give the output of sirius_param_files as an argument to this function, which contains a table of all the input and output paths for SIRIUS. Keep QC = TRUE, if QC files were used. If a suspect list is used, then SL = TRUE. Give path to the fragmentation tree directory. Candidates are the number of molecular formulas considered for further SIRIUS calculations. If SL = FALSE, keep the SL_path = NA
```R
run_sirius(files = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA_SiriusPandSL.csv', sep = ""), "./"), sep =""), 
               ppm_max = 5, 
               ppm_max_ms2 = 15, 
               QC = FALSE, 
               SL = TRUE, 
               SL_path = paste(input_dir, 'SL_Frag/', sep = ""),
               candidates = 30)
```
  
12. To run MetFrag, SIRIUS results are preprocessed to extract the best scored Molecular formula and its adduct. The Adduct annotation is important to run Metfrag. In this function, give the result directory and whether a suspect list was used or not.
```R
sirius_pproc <- sirius_postprocess(input_table[i, "ResultFileNames"], SL = TRUE)
```
                
13. Generate MetFrag parameter files. Give the SIRIUS processed results as input. Add the path of the MetFrag_AdductTypes.csv. Add the path of the suspect list InChIKeys generated by python function ```slist_metfrag``` for Metfrag.
```R
met_param <- metfrag_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/MS1DATAsirius.csv", sep = ""), "./"), sep =""), 
                               result_dir = input_table[i, "ResultFileNames"],
                               input_dir,
                               adducts = paste(input_dir, "MetFrag_AdductTypes.csv", sep = ""), 
                               sl_mtfrag = paste(input_dir, "SL_metfrag.txt", sep = ""), 
                               SL = TRUE,
                               ppm_max = 5, 
                               ppm_max_ms2= 15)
```
  
14. Run MetFrag.
```R
run_metfrag(met_param = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/metparam_list.txt", sep = ""), "./"), sep =""),
               MetFragjarFile = paste(input_dir, "MetFragCommandLine-2.4.8.jar", sep =""))
```

### Python
Once R workflow has been run for each file and there are results generated, python is used to preprocess that data. Follow the jupyter notebook Workflow_Python_Script.ipynb.
1. Import the python function file
```python
from Workflow_Python_Functions import *
```
2. Define input directory, keep all files in same directory and scripts so getwd works.
```python
input_dir = os.getcwd()+'/'
input_dir
```
3. Store suspect list in a variable to be used in the functions later.
```python
slistcsv = input_dir + "SkeletonemaSuspectListV1.csv"
```
4. This function states SIRIUS Post processing part 2. Here additionally maximum common substructures (MCSS) from the top scoring candidates are calculated by RDKit.
```python
sirius_postProc2(input_dir, 
                 input_tablecsv = input_dir + "input_table.csv")
```                 
5. This function performs metfrag post processing; from selecting the best candidate from KEGG and PubChem results, and also calcultaing MCSS for each feature.
```python
metfrag_postproc(input_dir, 
                 input_tablecsv = input_dir + "input_table.csv")
```
6. Combine results from different files. Note that source can be either all_insilico, SIRIUS, or MetFrag. The result will be combined from different files either with SIRIUS or MetFrag or both.
```python
combine_insilico(input_dir, 
                 input_tablecsv = input_dir + "input_table.csv",
                Source = "all_insilico")
```
7. Curation of results from SIRIUS
sirius_curation(input_dir, 
                 siriuscsv = input_dir + "MetabolomicsResults/SIRIUS_combined.csv", 
                 sl = True)
8. Curation of results from MetFrag
metfrag_curation(input_dir, 
                 metfragcsv = input_dir + "MetabolomicsResults/MetFrag_combined.csv", 
                 sl = True)

9. combineSM(input_dir, 
          metfragcsv = input_dir + 'MetabolomicsResults/metfrag_curated.csv', 
          siriuscsv = input_dir + 'MetabolomicsResults/sirius_curated.csv')
10. check each mzml file and each database csv result file; perform post processing
spec_postproc(input_dir, 
             Source = "all")

11. # combine all spectral databases for each mzml file
combine_specdb(input_dir)

12. # combine all spectral databases for all mzml file
combine_allspec(input_dir)

13. # only keep good scoring spectral database results
scoring_spec(input_dir, 
             spec_file = input_dir + 'MetabolomicsResults/SD_post_processed_combined_results.csv')

14. This function performs suspect list screening against spectal databases. The inputs can be db(str): all, gnps, mbank, hmdb, gm(gnps, mbank), hg(hmdb and gnps), hm(hmdb and mbank). It runs tanoimoto similarity score to between compounds from the results from spectral DBs and suspect list.
```python
suspectListScreening(input_dir, 
                     slistcsv, 
                     SpectralDB_Results = input_dir + 'MetabolomicsResults/scoredSpecDB.csv', 
                     db = "all")
```
15. combine_CuratedR prioritizes in the following manner: sirius>metfrag(kegg)>gnps>mbank>metfrag(pubchem)>hmdb + Suspectlist
```python
specDB_Curation(input_dir, 
                combinedx = input_dir + 'MetabolomicsResults/SpecDBvsSL.csv', 
                sl = True, 
                db = "all")

16. 
```python
combine_CuratedR(input_dir, 
                 combinedSDBs = input_dir + 'MetabolomicsResults/curatedSDB.csv', 
                 combinedSMs = input_dir + 'MetabolomicsResults/combinedSM.csv', 
                 data_type = "standards")
```
17. This function is to check the vailidity of SMILES of the tentative candidate.
```python
checkSMILES_validity(input_dir, 
                     resultcsv = input_dir + 'MetabolomicsResults/final_curation_without_classes.csv')
```
18. Classification function is performed for any feature which is not predicted by SIRIUS/CANOPUS. Such features have classifcation based on their SMILES, which are taken as input by ClassyFire to generate chemixl classes. Since both CANOPUS and ClassyFire use ChemONT, so the results are comparable.
```python
classification(input_dir, 
               resultcsv = input_dir + 'MetabolomicsResults/final_curation_with_validSMILES.csv')
```
## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)

