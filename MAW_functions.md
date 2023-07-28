# Tutorial of MAW-R Workflow according to the CWL file (Workflow_R_Script_all.r)

## Input files and Directories

An input directory (/input_dir) should have the following files.
1. One LCMS-2 spectra .mzML file
2. hmdb.rda (10.5281/zenodo.7519270)
3. mbankNIST.rda (10.5281/zenodo.7519270)
4. gnps.rda (10.5281/zenodo.7519270)
5. COCONUT.csv (10.5281/zenodo.7704937)

## MAW-R

Follow the R script: MAW/cwl/Workflow_R_Script_all.r. It starts with loading all the necessary libraries and then all the functions are written in the same file. At the end of the script, at line 2165, the actual executable script starts. The beginning refers to the predefining the parallelisation of the some functions. <br> The script is written in a way to run via Rscript command line, with the following arguments:

```R
#define arguments
args = commandArgs(trailingOnly = TRUE)

# input directory files
mzml_file <- args[1] # path to the mzml input file
gnps_file <- args[2] # path to the gnps.rda downloaded from Zenodo
hmdb_file <- args[3] # path to the hmdb.rda downloaded from Zenodo
mbank_file <- args[4] # path to the mbankNIST.rda downloaded from Zenodo
file_id <- args[5] # any given file ID
ppmx = as.numeric(args[6]) # ppm value for MS2 peak range
collision_info = as.logical(args[7]) # whether files have collision energy information
db_name = args[8] # any given name of the local database for MetFrag
db_path = args[9] # the path to the local database for MetFrag
```

Finally, generate a result directory using:
```R
mzml_result <- str_remove(basename(mzml_file), ".mzML")
dir.create(mzml_result)
```

### 1. MS2 data pre-processing
Now, we can start with the MS2 data pre-processing which includes reading the mzML files, removing empty spectra and extracting precursor m/z values. Here x is the mzML input file.
```R
# read mzML file and create output directory
spec_pr <- spec_Processing(mzml_file, mzml_result)
```
This results into two files. One file saves the pre-processed spectra in mzML format, accessed by spec_pr[[1]] (e.g: /opt/workdir/data/File1/processedSpectra.mzML) and other file keeps a list of precursor m/z(s) as a txt file, accessed by spec_pr[[2]] (e.g: /opt/workdir/data/File1/premz_list.txt)

### (Optional) 2. Download Spectral Databases (generally not required and is not part of the workflow)

MAW integrates three open-source spectral Databases: GNPS, HMDB and MassBank. Each of these databases can be downloaded using their specific file formats.
1. GNPS (ALL_GNPS) can be downloaded via [this link](https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.mgf) in mgf format.
2. MassBank_NIST can be downloaded with [this link](https://github.com/MassBank/MassBank-data/releases) in msp format.
3. HMDB (All Spectra Files)can be downloaded with this link(https://hmdb.ca/downloads) in xml format.
<br>
You can download the following versions of these databases from 10.5281/zenodo.7519270. <br>
GNPS saved at 2023-01-09 15:24:46<br>
HMDB saved at 2023-01-09 14:35:46 with release version Current Version (5.0)<br>
MassBank saved at 2022-01-09 15:10:52 with release version 2022.12 as mbankNIST.rda <br>
Store these .rda files in the same directory as your input mzml file e.g: opt/workdir/data

### 3. Spectral Database Dereplication

The module called ```spec_dereplication_file``` is used to perform Spectral database dereplication against the input data present in the mzML file. The function can be used in the following way:

```R
df_derep <- spec_dereplication_file(mzml_file = mzml_file,
                                    pre_tbl = paste(mzml_result, "/premz_list.txt", sep = ""),
                                    proc_mzml = paste(mzml_result, "/processedSpectra.mzML", sep = ""),
                                    db = "all",
                                    result_dir = mzml_result,
                                    file_id,
                                    no_of_candidates = 50,
                                    ppmx)
```
Here, mzml_file is the input file, pre_tbl is the list of precursor m/z(s) saved in a text file and proc_mzml is the pre-processed data saved as mzML file. db can be either "gnps", "hmdb", "mbank" or "all". no_of_candidates define how many candidates will be considered.

The function first performs further pre-processing steps needed for spectral database dereplication, such as removing low intensity peaks from MS2 fragmentation spectra, normalisation of the peak intensities and removal of any peaks higher or equal to precursor m/z value. This is performed for both the input spectra and the spectra in all the databases which matches the precursor m/z of the input mzML files. Once the preprocessing is done, the input spectra and the database candidate spectra are matched resulting in a score matrix.

Here, based on the score observation, GNPS >= 0.85, MassBank and HMDB >= 0.75 is the threshold for cosine similarity score. Another function within this module calculates the difference between matching peaks intensity and individual peak m/z. Based on the above functions, each precursor m/z from each file has a GNPS, MassBank and HMDB .csv files based on the selection of the database. e.g: /opt/workdir/data/File1/spectral_dereplication/MassBank/mbank_results_for_file_1M182RNAID1.csv would look like this

|no.| MBmax_similarity |MBmzScore|MBintScore|MQMatchingPeaks|MBTotalPeaks|mQueryTotalPeaks|MBFormula|MBSMILES|MBspectrumID|MBcompound|Source|
|---|--------------------|-----------|------------|---------------|--------------|----------------|----|----------|------|---|---|
| 1 | 0.971473979757729  |    0.933333333333333    |0.951921613619267|7|7|8|C9H11NO3|c1cc(ccc1C[C@@H](C(=O)O)N)O|MSBNK-UFZ-UA005601| L-Tyrosine|MassBank|

Here the name of the file indicates the file id which is file_1, the precursor m/z which is M182 and the retention time is not goven hence RNA, with the id ID1. The .csv file contains information on MBmax_similarity which is the matching spectra cosine similarity score, MBmzScore which is the score for matching m/z of the peaks and ranges from 0-1, MBintScore which is the score for matching intensity of the peaks and ranges from 0-1, MQMatchingPeaks which is the number of matching peaks between MassBank and Query spectra, and other columns names are self-explanatory. This form of result follows for GNPS and HMDB as well.

### 4. MAW-R Compound Databases Dereplication

After the pre-processing steps from section 1, we can also do the compound database dereplication. This module is also divided into further independent functions. Here we use MetFrag as a tool for insilico fragmentation spectral matching, but the workflow also generates SIRIUS input .ms files. <br>
Metfrag .txt parameter file
```
PeakListPath = File1/insilico/peakfiles_ms2/Peaks_01.txt
IonizedPrecursorMass = 152.995452880859
PrecursorIonMode = -1
IsPositiveIonMode = False
MetFragDatabaseType = LocalCSV
LocalDatabasePath = COCONUT_Jan2022.csv
DatabaseSearchRelativeMassDeviation = 5
FragmentPeakMatchAbsoluteMassDeviation = 0.001
FragmentPeakMatchRelativeMassDeviation = 15
MetFragCandidateWriter = CSV
SampleName = 1_id_endo_negM153R40ID1_mz_152.995452880859_rt_40.0539309
ResultsPath = File1/insilico/MetFrag/coconut/
MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter
MetFragPostProcessingCandidateFilter = InChIKeyFilter
MaximumTreeDepth = 2
NumberThreads = 1
```

SIRIUS .ms file
```
>compound file_1M182RNAID1 ## ID
>parentmass 151.03643799   ## precursor m/z
>charge -1                 ## charge, pos or neg
>rt NAs              ## median retention time in seconds
>collision 15eV            ## ms2 simply can be written here instead of collision energy, this is the fragmentation peak list
65.1203 26499.103516
67.302422 30616.269531
71.056732 26360.367188
...
136.075729 7684037
147.043976 1080712.25
165.054611 3833276.5
182.081146 575214.9375  ### first column is m/z and second column is intensity
```

In order to create such a file, further pre-processing and data extraction is required from mzML format to be stored in txt or .ms. To do that, ```ms2_peaks``` function is used, which extracts ms2 fragmentation peaks for all of the precursor m/z. If there are more MS2 spectra associated with one precursor m/z, the function combines these MS2 spectra and extracts all fragmentation peaks. These peaks are stored as txt (e.g: "/usr/MAW/Tutorial/File1/insilico/peakfiles_ms2/Peaks1.txt")
```R
# Extract MS2 peaks
spec_pr2 <- ms2_peaks(pre_tbl = paste(mzml_result, "/premz_list.txt", sep = ""),
                      proc_mzml = paste(mzml_result, "/processedSpectra.mzML", sep = ""),
                      result_dir = mzml_result,
                      file_id)
# Extract MS1 peaks
ms1p <- ms1_peaks(x = paste(mzml_result,'/insilico/MS2DATA.csv', sep = ""),
                    y = NA, 
                    result_dir = mzml_result,
                    QCfile = FALSE)
```
The result from this function is not only the MS2 Fragmentation peak lists, but also a table of precursor m/z along with other elements of each feature. e.g: the file "/opt/workdir/data/File1/insilico/MS2DATA.csv" contains:

|no.|       id_X     |   premz    |  rtmed | rtmean |int|col_eng|pol|ms2Peaks|ms1Peaks|
|---|----------------|------------|--------|--------|---|-------|---|-----------|-----------|
| 1 |file_1M182RNAID1|182.081|NA|NA| NA|   NA |neg|./Example_Tyrosine/insilico/peakfiles_ms2/Peaks_01.txt|NA|

Here, id represents a unique id for individual feature (precursor m/z and rt median). This ID corresponds with the one generated for spectral database dereplication. Last column represents the path of the MS2 fragmentation peaks which can later be used to retrieved the ms2 peaks.<br>

(OPTIONAL) Now, to create the ms files, we have extracted all relevant information. To write the .ms files, function called ```sirius_param``` is used. MS1DATA.csv has all the information needed to write the .ms file.
```R
#prepare sirius parameter files
sirius_param_files <- sirius_param(x = paste(mzml_result,'/insilico/MS1DATA.csv', sep = ""),
                       result_dir = mzml_result,
                       SL = FALSE, 
                       collision_info)
```
The resulting table looks like this:
|no.|sirius_param_file|   outputNames    | isotopes |
|---|----------------|------------|------------|
| 1 |/opt/workdir/data/File1/insilico/SIRIUS/1_NA_iso_MS1p_182_SIRIUS_param.ms|/opt/workdir/data/File1/insilico/SIRIUS/1_NA_iso_MS1p_182_SIRIUS_param.json|NA|

To generate txt files for metfrag following command is used:
```R
metfrag_param(x= paste(mzml_result,'/insilico/MS1DATA.csv', sep = ""), 
                result_dir = mzml_result, 
                db_name,
                db_path, 
                ppm_max = 5, 
                ppm_max_ms2= 15)
```
After that a JSON file is created with all the outputs and their relevant paths present in the JSON file.

# Tutorial on usage of MetFrag within CWL

CWL integrates the docker image for MetFrag to run MetFrag using the parameters in JSON file generated by MAW-R.

```yaml
cwlVersion: v1.0
class: CommandLineTool

baseCommand: ['java', '-jar', '/usr/src/myapp/MetFragCommandLine-2.5.0.jar']

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/maw-metfrag_2.5.0:1.0.5
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - entryname: metfrag.inputs
      entry: |-
        PeakListPath = $(inputs.PeakList.path)
        IonizedPrecursorMass =  $(inputs.IonizedPrecursorMass)
        PrecursorIonMode = $(inputs.PrecursorIonMode)
        MetFragDatabaseType = LocalCSV
        LocalDatabasePath = $(inputs.LocalDatabasePath.path)
        DatabaseSearchRelativeMassDeviation = 5
        FragmentPeakMatchAbsoluteMassDeviation = 0.001
        FragmentPeakMatchRelativeMassDeviation = 15
        MetFragCandidateWriter = CSV
        SampleName = $(inputs.SampleName)
        #ResultsPath = $(runtime.outdir)/metfrag
        ResultsPath = $(runtime.outdir)
        MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter
        MetFragPostProcessingCandidateFilter = InChIKeyFilter
        MaximumTreeDepth = 2
        NumberThreads = 1


inputs: # additional inputs for all files; make them to show certain paths
  PeakList:
      type: File
  IonizedPrecursorMass:
      type: string
  PrecursorIonMode:
      type: int
  LocalDatabasePath:
      type: File
  SampleName:
      type: string

arguments: 
  - $(runtime.outdir)/metfrag.inputs

outputs:
  metfrag_candidate_list:
    type: File
    outputBinding:
        glob: "$(runtime.outdir)/*.csv"
```

The JSON file fills in the metfrag.inputs for each parameter file, and only few of these parameters are dynamic which are then listed in inputs section:   PeakList which is MS2 peak list, IonizedPrecursorMass which precursor mass, PrecursorIonMode which is the ionisation mode, LocalDatabasePath which is path to COCONUT database, and SampleName which self given name. <br>

Alternatively you can run the Run_MetFrag.r instead of CWL, but you will need the table from ```metfrag_param``` and the parameter files.

# Tutorial of MAW-Python Workflow with CWL

MAW-Python is dependent on the results obtained from MAW-R and MAW-MetFrag. Following the Workflow_Python_Script_all.py, the beginning calls all the packages, then the functions are defined. Around the line 3002, the arguments for the python3 command are written.

```python
# Define the command-line arguments
parser = argparse.ArgumentParser(description='MAW-Py')
parser.add_argument('--file_id', type=str, help='file_id')
parser.add_argument('--msp_file', type=str, help='path to spec result CSV file')
parser.add_argument('--gnps_dir', type=str, help='path to GNPS directory')
parser.add_argument('--hmdb_dir', type=str, help='path to HMDB directory')
parser.add_argument('--mbank_dir', type=str, help='path to MassBank directory')
#name it metfrag_candidate
parser.add_argument('--metfrag_candidate_list', type=str, action='append', help='path to MetFrag candidate table CSV file')
parser.add_argument('--ms1data', type=str, help='path to MS1 data CSV file')
parser.add_argument('--score_thresh', type=float, default=0.75, help='score threshold for MetFrag results (default: 0.75)')
# Parse the command-line arguments
args = parser.parse_args()
```

### 1.  MAW-R Results Post-processing
The first function post processes the results from GNPS, MassBank and HMDB and writes new files named liked proc.csv. For GNPS, the Compound names and SMILES need to be curated, so MAW-python uses RDKit and PubChemPy to retrieve this information. For HMDB, the original results only have the ID, so a sdf file is downloaded from this [link](https://hmdb.ca/downloads) to extract further information. It also removes any low scoring candidates.
```python
print("spec_postproc starts")
msp_file_df = spec_postproc(msp_file, gnps_dir, hmdb_dir, mbank_dir, file_id)
```
### 2  MetFrag Results Post-processing
The second step is to post-process the results from MetFrag, where any low scroing candidates are simply removed.
```python
print("metfrag_postproc starts")
ms1data_df = metfrag_postproc(ms1data, metfrag_candidate_list, file_id, score_thresh)
```
### 3.  Candidate Selection
Candidate selection function creates a folder /CandidateSelection/ which stores a file called ChemMN.tsv for each feature. This .tsv file can be imported as a network in Cytoscape to visualize a chemical similarity network. Keep Source and Target nodes as the Names of the candidates such as M_1(MassBank candidate ranked 1) and use tanimoto score as edge. The second type of file is a ranked list of candidates generated with MAW called sorted_candidate_list.csv. If there are more files, merge the results from all mergedResults-with-one-Candidates.csv files from each file.
```python
print("CandidateSelection starts")
CandidateSelection_SimilarityandIdentity_Metfrag(file_id = file_id, msp_file = msp_file_df, 
ms1data = ms1data_df, standards = False)
```
### 4. Classification

Classification function is performed for all features. Such features have classification based on their SMILES, which are taken as input by ClassyFire to generate chemical classes. 
```python
classification(resultcsv = file_id + "_mergedResults-with-one-Candidates.csv")
```