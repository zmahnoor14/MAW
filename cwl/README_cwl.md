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
5. hmdb_dframe_str.csv

## Run maw.cwl to collect provenance:

To run the workflow without collecting provenance is:
```
cwltool --relax-path-check --cachedir cache --outdir run12  maw.cwl --r_script Workflow_R_Script_all.r --python_script Workflow_Python_Script_all.py --mzml_files ./mzml_Files 2> run12.log
```

Make sure you have the paths entered correctly. Your working directory where you have kept maw.cwl, should also have both scripts and the maw-r.cwl and maw-py.cwl. In your working directory you can create another directory called mzml_Files folder which shoul have all the input files mention in the section above.

To run the workflow **with** provenance:
```
cwltool --relax-path-check --cachedir cache --outdir run12 --provenance /any/writeable/directory/to/store/provenance maw.cwl --r_script Workflow_R_Script_all.r --python_script Workflow_Python_Script_all.py --mzml_files ./mzml_Files 2> run12.log
```
