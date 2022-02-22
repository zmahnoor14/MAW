#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load Libraries
library(Spectra)
library(MsBackendMgf)
library(MsBackendHmdb)
library(MsCoreUtils)
library(MsBackendMsp)
library(readr)
library(dplyr)
# 3 dependencies for latest MassBank version
library(rvest)
library(stringr)
library(xml2)

# ---------- Script ----------
# input directory
input_dir <- paste(getwd(), "/", sep = '')

# load the functions file
source(file = paste(input_dir, "Workflow_R_Functions.r", sep = ''))

# load db spectra objects [gnps, hmdb, mbank]
#download_specDB(input_dir, db = "all")

# OR load without using function 
load(file = paste(input_dir,"gnps.rda", sep = ""))
load(file = paste(input_dir,"hmdb.rda", sep = ""))
load(file = paste(input_dir,"mbank.rda", sep = ""))

# create result folders and store in a dataframe
input_table <- data.frame(ms2_rfilename(input_dir))

#cam_funcMode(path = paste(input_dir, "QC", sep =""), pattern = "common")
#merge_qc(path = paste(input_dir, "QC", sep =""))
#cam_func(path = "QC/", f = "DS200212_Scost_QC_280k_pos.mzML", mode = "pos")
#cam_func(path = "QC/", f = "DS200212_Scost_QC_280k_neg.mzML", mode = "neg")
for (i in 1:nrow(input_table)){
    # if a certain phrase is present in the data files e.g: pos, then take the pos CAMERA
    if (grepl("SC_full_PRM_pos", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_pos.csv"
    }
    if (grepl("SC_full_PRM_neg", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_neg.csv"
    }
    if (grepl("DS200309_Scost_QC_70k_pos_PRM.mzML", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/posCAMERAResults_DS200212_Scost_QC_280k_pos.csv"
    }
    if (grepl("DS200309_Scost_QC_70k_neg_PRM.mzML", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/negCAMERAResults_DS200212_Scost_QC_280k_neg.csv"
    }
}

for (i in 1:nrow(input_table)){
    # Preprocess and Read the mzMLfiles
    spec_pr <- spec_Processing(as.character(input_table[i, "mzml_files"]))
    # Extract spectra
    sps_all <- spec_pr[[1]]
    # Extract precursor m/z
    pre_mz<- spec_pr[[2]]
    
    for (a in pre_mz){
        #perform dereplication with all dbs
        df_derep <- spec_dereplication(a, db = "all", 
                                       result_dir = input_table[i, "ResultFileNames"],
                                       file_id = input_table[i, "File_id"], input_dir, 
                                       ppmx = 15)
    }
    # Extract MS2 peak lists
    spec_pr2 <- ms2_peaks(spec_pr, input_table[i, "ResultFileNames"])
    
    # Extract MS1 peaks or isotopic peaks
    ms1p <- ms1_peaks(spec_pr2, input_table[i, "qcCAM_csv"], input_table[i, "ResultFileNames"], QC = TRUE)
    
    #prepare sirius parameter files
    sirius_param_files <- sirius_param(ms1p, 
                                       result_dir = input_table[1, 'ResultFileNames'], 
                                       SL = TRUE)
    
    # Run sirius
    run_sirius(files=sirius_param_files , 
               ppm_max = 5, 
               ppm_max_ms2 = 15, 
               QC = TRUE, 
               SL = TRUE, 
               SL_path = paste(input_dir, 'ScostSLS/', sep = ""),
               candidates = 30)
    
    # Post process Sirius results and extract adducts for MetFrag
    sirius_pproc <- sirius_postprocess(input_table[i, "ResultFileNames"], SL = TRUE)
    
    # prepare Metfrag parameter files
    met_param <- metfrag_param(sirius_pproc, 
                               result_dir = input_table[1, "ResultFileNames"],
                               input_dir, 
                               adducts = paste(input_dir, "MetFrag_AdductTypes.csv", sep = ""), 
                               sl_mtfrag = paste(input_dir, "sl_metfrag.txt", sep = ""), 
                               SL = TRUE)
    # run metfrag
    for (files in met_param){
        system(paste("java -jar",  paste(input_dir, "MetFrag2.4.5-CL.jar", sep = ''), files))
    }
}


