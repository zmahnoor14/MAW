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
input_dir

# load the functions file
source(file = paste(input_dir, "Worlfow_R_Functions.r", sep = ''))

# load db spectra objects [gnps, hmdb, mbank]
#download_specDB(input_dir, db = "all")
#load_specDB(input_dir, db = "gnps") # doesnt work



# OR load without using function 
load(file = paste(input_dir,"gnps.rda", sep = ""))
load(file = paste(input_dir,"hmdb.rda", sep = ""))
load(file = paste(input_dir,"mbank.rda", sep = ""))

# create result folders and store in a dataframe
input_table <- data.frame(ms2_rfilename(input_dir))
input_table

spec_pr <- spec_Processing(as.character(input_table[1, "mzml_files"]))
sps_all <- spec_pr[[1]]
pre_mz<- spec_pr[[2]]
#x <- as.numeric(pre_mz)

df_derep <- spec_dereplication(pre_mz[1], 
                                       db = "all", 
                                       result_dir = input_table[1, "ResultFileNames"],
                                       file_id = input_table[1, "File_id"],
                                       input_dir,
                                       ppmx = 15)

spec_pr2 <- ms2_peaks(spec_pr, input_table[1, "ResultFileNames"])

cam_funcMode(path = paste(input_dir, "QC", sep =""), pattern = "common")

merge_qc(path = paste(input_dir, "QC", sep =""))

cam_func(path = "/Users/mahnoorzulfiqar/Project_Scost/QC/", f = "DS200212_Scost_QC_280k_pos.mzML", mode = "pos")
cam_func(path = "/Users/mahnoorzulfiqar/Project_Scost/QC/", f = "DS200212_Scost_QC_280k_neg.mzML", mode = "neg")

for (i in 1:nrow(input_table)){
    # if a certain phrase is present in the data files e.g: pos, then take the pos CAMERA
    if (grepl("SC_full_PRM_pos", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_pos.csv"
    }
    if (grepl("SC_full_PRM_neg", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_neg.csv"
    }
    if (grepl("DS200309_Scost_QC_70k_pos_PRM.mzML", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/posCAMERA_Results_DS200212_Scost_QC.csv"
    }
    if (grepl("DS200309_Scost_QC_70k_neg_PRM.mzML", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/negCAMERA_Results_DS200212_Scost_QC.csv"
    }
}
input_table

ms1p <- ms1_peaks(x = spec_pr2, 
                  y = input_table[1, 'qcCAM_csv'], 
                  result_dir = input_table[1, 'ResultFileNames'], 
                  QC = TRUE)

sirius_param_files <- sirius_param(ms1p, 
                                   result_dir = input_table[1, 'ResultFileNames'], 
                                   SL = TRUE)

sirius_param_files



adducts <- paste(input_dir, "MetFrag_AdductTypes.csv", sep = '')

# ---------- multiple files ----------

for (i in 1:nrow(input_table)){
    spec_pr <- spec_Processing(as.character(input_table[i, "mzml_files"]))
    sps_all <- spec_pr[[1]]
    pre_mz<- spec_pr[[2]]
    
    for (a in pre_mz){
        
        df_derep <- spec_dereplication(a, 
                                       db = "all", 
                                       result_dir = input_table[i, "ResultFileNames"],
                                       file_id = input_table[i, "File_id"],
                                       input_dir,
                                       ppmx = 15)
    }
    
    spec_pr2 <- ms2_peaks(spec_pr, input_table[i, "ResultFileNames"])
    spec_pr3 <- ms1_peaks(spec_pr2, y = NA, result_dir = input_table[i, "ResultFileNames"], QCfile = FALSE)
    sirius_pfiles <- sirius_param(spec_pr3, result_dir = input_table[i, "ResultFileNames"], SL = FALSE)
    
    for (b in nrow(sirius_pfiles)){
        system(paste("sirius --input", sirius_pfiles[b, "sirius_param_file"], 
                     "--output", sirius_pfiles[b, "outputNames"],
                    "formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates 30 --ppm-max 5 --ppm-max-ms2 15 structure --database ALL canopus", 
                     sep = " "))
    }
    
    sirius_pproc <- sirius_postprocess(input_table[i, "ResultFileNames"], SL = FALSE)
    
    met_param <- metfrag_param(sirius_pproc, result_dir = input_table[i, "ResultFileNames"],
                               input_dir, adducts, sl_mtfrag = NA, SL = FALSE)
    for (files in met_param){
        system(paste("java -jar",  paste(input_dir, "MetFrag2.4.5-CL.jar", sep = ''), files))
    }
}






