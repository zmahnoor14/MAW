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
options(warn=-1)

# Track Time 
start_time <- Sys.time()

# ---------- Script ----------
# input directory
input_dir <- paste(getwd(), "/", sep = '')
input_dir

#input_dir <- "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/"
#input_dir

# load the functions file
source(file = paste(input_dir, "Workflow_R_Functions.r", sep = ''))

# load the functions file
# source(file = '/Users/mahnoorzulfiqar/OneDriveUNI/MAW/Workflow_R_Functions.r')

# downloading spectral libraries; do NOT run
# load db spectra objects [gnps, hmdb, mbank]
# download_specDB(input_dir, db = "all")

# OR load the database rda objects 
load(file = paste(input_dir,"gnps.rda", sep = ""))
load(file = paste(input_dir,"hmdb.rda", sep = ""))
load(file = paste(input_dir,"mbank.rda", sep = ""))

#gnps

#mbank

#hmdb

# Run the first function; this creates a dataframe of your input files, their result directories 
# and gives an id to each input file; stores the table in directory as a csv filr
input_table <- data.frame(ms2_rfilename(input_dir))
input_table

# Do not run these functions 
#cam_funcMode(path = paste(input_dir, "QC", sep =""), pattern = "common")
#merge_qc(path = paste(input_dir, "QC", sep =""))
#cam_func(path = "QC/", f = "DS200212_Scost_QC_280k_pos.mzML", mode = "pos")
#cam_func(path = "QC/", f = "DS200212_Scost_QC_280k_neg.mzML", mode = "neg")

## Here added are the QC files, run the function as it is

# add QC files to respective pos and neg file
for (i in 1:nrow(input_table)){
    # if a certain phrase is present in the data files e.g: pos, then take the pos CAMERA
    if (grepl("SC_full_PRM_pos", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_pos.csv"
    }
    if (grepl("SC_full_PRM_neg", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/Combined_Camera_neg.csv"
    }
    if (grepl("DS200309_Scost_QC_70k_pos_PRM", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/posCAMERAResults_DS200212_Scost_QC_280k_pos.csv"
    }
    if (grepl("DS200309_Scost_QC_70k_neg_PRM", input_table[i, "mzml_files"], fixed=TRUE)){
        input_table[i, "qcCAM_csv"] <- "./QC/negCAMERAResults_DS200212_Scost_QC_280k_neg.csv"
    }
}

input_table

#Run a loop for the two files



for (i in 1:nrow(input_table)){
    
    
    # Preprocess and Read the mzMLfiles
    spec_pr <- spec_Processing(as.character(input_table[i, "mzml_files"]), input_table[i, "ResultFileNames"])
    
    
    # Extract spectra
    sps_all <- spec_pr[[1]]
    # Extract precursor m/z
    pre_mz<- spec_pr[[2]]

    #perform dereplication with all dbs
    df_derep <- spec_dereplication(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "./"), sep =""), 
                                   proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "./"), sep =""),
                                   db = "all", 
                                   result_dir = input_table[i, "ResultFileNames"],
                                   file_id = input_table[i, "File_id"], 
                                   input_dir, 
                                   ppmx = 15)
    
    
    # Extract MS2 peak lists
    spec_pr2 <- ms2_peaks(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "./"), sep =""), 
                          proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "./"), sep =""),
                          input_dir,
                          input_table[i, "ResultFileNames"],
                         file_id = input_table[i, "File_id"]) 
    
    # Extract MS1 peaks or isotopic peaks
    ms1p <- ms1_peaks(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS2DATA.csv', sep = ""), "./"), sep =""), 
                      y = input_table[i, "qcCAM_csv"], 
                      input_table[i, "ResultFileNames"], 
                      input_dir, 
                      QC = TRUE)
    
    #prepare sirius parameter files
    sirius_param_files <- sirius_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA.csv', sep = ""), "./"), sep =""), 
                                       result_dir = input_table[i, 'ResultFileNames'], 
                                       input_dir,
                                       SL = TRUE)
    
    # Run sirius
    run_sirius(files = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA_SiriusPandSL.csv', sep = ""), "./"), sep =""), 
               ppm_max = 5, 
               ppm_max_ms2 = 15, 
               QC = TRUE, 
               SL = TRUE, 
               SL_path = paste(input_dir, 'ScostSLS/', sep = ""),
               candidates = 30)
    
    
    
    # Post process Sirius results and extract adducts for MetFrag
    sirius_pproc <- sirius_postprocess(input_table[i, "ResultFileNames"], SL = TRUE)
    
    
    
    # prepare Metfrag parameter files
    met_param <- metfrag_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/MS1DATAsirius.csv", sep = ""), "./"), sep =""), 
                               result_dir = input_table[i, "ResultFileNames"],
                               input_dir,
                               adducts = paste(input_dir, "MetFrag_AdductTypes.csv", sep = ""), 
                               sl_mtfrag = paste(input_dir, "sl_metfrag.txt", sep = ""), 
                               SL = TRUE,
                               ppm_max = 5, 
                               ppm_max_ms2= 15)
    
    
    # run metfrag
    run_metfrag(met_param = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/metparam_list.txt", sep = ""), "./"), sep =""),
                input_dir)
    
    
}

end_time <- Sys.time()
print(end_time - start_time)


