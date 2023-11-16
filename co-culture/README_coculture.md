# Co-culture workflow
This workflow can be used to pre-process and analyse MS1 data stemming from a two species, untargeted, comparative metabolomics liquid-chromatography tandem mass spectrometry approach. It takes .mzML files as input and gives processed feature lists and statistical analyses of the MS level 1 data.

The workflow was developed using a co-culture experiment between the two microalgal species _S. marinoi_ and _P. parvum_

## HOW TO:
1. The _MS1 workflow.R_ file uses the functions defined in the _MS1 workflow functions.R_ file to run the workflow with the following pre-processing steps using XCMS functionalities:
	- data_preparation: the .mzML files are loaded and an XCMSnExp object is created for the MS data
	- chromatogram_qc: chromatograms are plotted for preliminary data inspection
	- peak_detection: chromatographic peaks are being detected 
	- grouping_1: detected peaks are being group throughout the samples
	- rt_correction: retention time correction is performed and plotted 
	- grouping_2: peaks are being grouped again after retention time correction
	- feature_extraction: the feature infomation is extracted from the XCMSnExp object
	- feature_transformation: a feature list is created with transformed data and imputed missing values
	- bina_list_creation: a binary list is created for subsequent statistical analysis

2. The inputs to be defined when using in CL mode are as follows:
	- files: list of .mzMl files
	- phenodata: .csv file of phenodata, including (1) sample_name, (2) sample_group, (3) sample_description
	- result_dir_name: name for directory where the results are stored
	- plots_dir_name: name for directory where the plots are stored
	- chrom_run_sec: length of chromatography run in seconds
	- msLevel: (optional) defaults to MS level 1 
	- CentWaveParam: parameters for peak picking in form of a CentWaveParam class from xcms. Default settings are given but optimization for data at hand is highly recommended
	- PeakDensityParam_gr: parameters for grouping in form of a PeakDensityParam class from xcms. Default settings are given as minFraction = 0.7, bw = 2.5
	- PeakGroupsParam_rt: parameters for retention time correction in form of a PeakGroupsParam class from xcms. Default settings are given as minFraction = 0.7, smooth = "loess", span = 0.5, family = "gaussian"
	- intensity_cutoff: for the creation of the binary list as absence/presence measure, a logarithmically transformed intensity cutoff, can be estimated from a histogram of the feature table

3. For in-depth statistical analysis, the _/Statistics/_ directory contains RScripts for several statistical methods. The objects used are the previously produced .csv tables of feat_list and bina_list. The negative and positive feature tables are combined and analysis is performed within the conditions (EXO and ENDO) combined and on species level. Statistical methods used include:
	- PCA
	- Diversity measures (Shannon diversity index, ...)
	- PLS 
	- t-SNE
	- Variation partitioning 

4. By following the _linking.R_ script, MS1 precursor can be linked to MS2 fragementation data. The MS2 fragmentation can be annotated using [MAW](https://github.com/zmahnoor14/MAW)



### Additional files and info
some statistical function are from the [iESTIMATE](https://github.com/ipb-halle/iESTIMATE/blob/main/use-cases/radula-hormones/peak_detection_neg.r) repository


