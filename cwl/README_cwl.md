# Generate CWLProv RO Crate Bundle 

## CWL Files needed to run maw.cwl:
1. maw-r.cwl
2. maw-py.cwl

## Script files needed to run maw.cwl:
1. Workflow_R_Script_all.r
2. Workflow_Python_Script_all.py

## Input data files:

1. Example_Tyrosine.mzML

These files can be found on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6528931.svg)](https://doi.org/10.5281/zenodo.6528931)

2. gnps.rda
3. hmdb.rda
4. mbankNIST.rda
5. hmdb_dframe_str.csv # for now optional

## Run maw.cwl to collect provenance:

#### (Optional if you are using M1 processor/(linux/arm64))

Setup the linux/amd64 processor:
```shell
$ export DOCKER_DEFAULT_PLATFORM=linux/amd64
```
#### (for intel, directly run the following)
To run the workflow without collecting provenance is:
```shell
$ cwltool --relax-path-check --cachedir cache --outdir ./run_maw_rpy maw.cwl --r_script Workflow_R_Script_all.r --python_script Workflow_Python_Script_all.py --mzml_files mzml_Files/Example_Tyrosine.mzML --gnps_rda mzml_Files/gnps.rda --hmdb_rda mzml_Files/hmdb.rda --mbank_rda mzml_Files/mbankNIST.rda 2> run_maw_rpy.log
```

Make sure you have the paths entered correctly. Your working directory where you have kept maw.cwl, should also have both scripts and the maw-r.cwl and maw-py.cwl. In your working directory you can create another directory called mzml_Files folder which should have all the input files mentioned in the section above.

To run the workflow **with** provenance:
```shell
$ cwltool --relax-path-check --cachedir cache --outdir ./run_maw_rpy --provenance /any/writeable/directory/to/store/provenance maw.cwl --r_script Workflow_R_Script_all.r --python_script Workflow_Python_Script_all.py --mzml_files mzml_Files/Example_Tyrosine.mzML --gnps_rda mzml_Files/gnps.rda --hmdb_rda mzml_Files/hmdb.rda --mbank_rda mzml_Files/mbankNIST.rda 2> run_maw_rpy.log
```
The addition **with** provenance is to add --provenance /any/writeable/directory/to/store/provenance before the cwl file.
