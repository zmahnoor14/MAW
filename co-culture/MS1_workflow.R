###---- library ----
library(RColorBrewer)           # For colors
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(vegan)                  # For shannon diversity
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(circlize)               # For sunburst plot
library(plotrix)                # For sunburst plot
library(caret)                  # Swiss-army knife for statistics
library(pROC)                   # Evaluation metrics
library(PRROC)                  # Evaluation metrics
library(multiROC)               # Evaluation metrics
library(plotly)                 # For creating html plots
library(htmlwidgets)            # For creating html plots
library(stringr)              

### MS1_workflow_functions
# this is an examplary run of the workflow using the functions from the functions file
# if you want to run the workflow on your own data edit the input parameters for the functions accordingly

# edit for your local path to the functions file
source("C:/Users/abela/Documents/GitHub/CoCultureSmPp/MS1_workflow_XCMS/MS1_workflow_functions.R")

# ---- analysis for ENDO neg ----
# defining variable for functions to use
# edit files path to your local path to the files
raw_data_MS1_ENDO_neg <- list.files("C:/Users/abela/Documents/Uni_Jena/Masterarbeit/MAW-Co-culture/MS1_pos_neg", pattern = "neg_ENDO", full.names = TRUE)
phenodata <- read.csv("phenodata_ENDO_neg.csv")
result_dir_name <- "ENDO_neg_testing_results"
plots_dir_name <- "ENDO_neg_testing_plots"
ms1_params_ENDO_neg <- CentWaveParam(ppm=25, mzCenterFun="wMean", peakwidth=c(14, 59), 
                                     prefilter=c(3, 140), mzdiff=0.0155, snthresh=7, noise=0, 
                                     integrate=1, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                     fitgauss=FALSE, roiList=list(), roiScales=numeric())

# prepare the data 
msd <- data_preparation(files = raw_data_MS1_ENDO_neg, phenodata = phenodata, result_dir_name = "ENDO_neg_testing_results", 
                        plots_dir_name = "ENDO_neg_testing_plots", chrom_run_sec = 700, msLevel = NULL)
# plot chromatogram
chrom_msd <- chromatogram_qc(msd, phenodata, result_dir_name, plots_dir_name)

# perform peak detection on chromatographic peaks
ms1_data <- peak_detection(msd, CentWaveParam = ms1_params_ENDO_neg, result_dir_name, plots_dir_name)

# grouping of the peaks
ms1_data <- grouping_1(ms1_data, PeakDensityParam_gr = NULL, result_dir_name)

# perform retention time correction
ms1_data <- rt_correction(ms1_data, plots_dir_name, result_dir_name, PeakGroupsParam_rt = NULL, chromas_msd = chrom_msd)

# second grouping of the peaks
ms1_data <- grouping_2(ms1_data, PeakDensityParam_gr = NULL, result_dir_name)

# create a feature list 
feature_list <- feature_extraction(ms1_data, result_dir_name)

# prepare the feature list for further analyses
feature_list <- feature_transformation(feature_list, result_dir_name, plots_dir_name)

# create a binary absence/presence list for statistical analysis
bina_list <- bina_list_creation(ms1_matrix = NULL, intensity_cutoff = 14, result_dir_name)

