# Installation on all systems

Please follow the instructions given on the MAW [README.md](https://github.com/zmahnoor14/MAW/tree/main?tab=readme-ov-file#installation) file. 

## Input files:
1. [cwl/maw-r.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-r.cwl)
2. [cwl/maw-py.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-py.cwl)
3. [cwl/maw-metfrag.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-metfrag.cwl)
4. [cwl/maw_single.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw_single.cwl)
5. [cwl/maw.cwl](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw.cwl)
6. [cwl/maw-inputs-multiple_ara.yaml](https://github.com/zmahnoor14/MAW/blob/main/cwl/maw-inputs-mutiple_ara.yaml)

The input file should contain paths to the following files:
1. [COCONUT Database](https://zenodo.org/record/7704937)
2. [Spectral Databases](https://zenodo.org/record/7519270): give the paths for each database separately in the yaml input file.
   a. [gnps_file](https://zenodo.org/records/7519270#:~:text=Download%20all-,gnps.rda,-md5%3Af5139892bdf216b54a6af9f1907f09ca)
   b. [hmdb_file](https://zenodo.org/records/7519270#:~:text=Download-,hmdb.rda,-md5%3Aa9dd9c1c3c023339a9e5c6fc4d4288e5)
   c. [mbank_file](https://zenodo.org/records/7519270#:~:text=Download-,mbankNIST.rda,-md5%3Af946d74093e819ca2078a7ef15a310c6)
4. [cwl/Workflow_R_Script_all.r](https://github.com/zmahnoor14/MAW/blob/main/cwl/Workflow_R_Script_all.r)
5. [cwl/Workflow_Python_Script_all.py](https://github.com/zmahnoor14/MAW/blob/main/cwl/Workflow_Python_Script_all.py)
6. Example files: I used the following: [VN_211016_Sc_st_PRM_neg.mzML](https://zenodo.org/records/7106205#:~:text=VN_211016_Sc_st_PRM_neg.mzML) and [VN_211016_Sc_QC_PRM_neg.mzML](https://zenodo.org/records/7106205#:~:text=VN_211016_Sc_QC_PRM_neg.mzML)


## MAC system

### System Specifications
* Chip: Apple M1 Max
* Memory: 32 GB

## Execute MAW on Mac Terminal
Make sure dopcker is running
```shell
cd cwl_input_dir
export DOCKER_DEFAULT_PLATFORM=linux/amd64
cwltool --parallel --outdir mac_outputs_parallel maw.cwl maw-inputs-mutiple_ara.yaml
```
## HPC Cluster 

### System Specifications
* Operating system: CentOS Linux release 7.6
* Processors: 140 nodes, each with 24 CPU cores

### Extra steps to installation on an HPC cluster
Load necessary modules: `module load tools/singularity`
Install nodejs and add CWL_SINGULARITY_CACHE  to ENV variable

## Execute MAW on Terminal
```shell
cd cwl_input_dir
cwltool --singularity --parallel --outdir ara_outputs maw.cwl maw-inputs-mutiple_ara.yaml
```
