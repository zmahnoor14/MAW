start.time <- Sys.time()

# ---------- Script ----------
# input directory
input_dir <- paste(getwd(), "/data", sep = "")
input_dir

# load the functions file
source(paste(getwd(), "/Workflow_R_Functions.r", sep = ""))

input_table <- data.frame(ms2_rfilename(input_dir))
input_table

for (i in 1:nrow(input_table)){
    
    #Preprocess and Read the mzMLfiles
    spec_pr <- spec_Processing(input_dir,
                               input_table[i, "mzml_files"], 
                               input_table[i, "ResultFileNames"])
    
    #perform dereplication with all dbs
    df_derep <- spec_dereplication_file(mzml_file = input_table[i, "mzml_files"],
                                   pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "."), sep =""), 
                                   proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "."), sep =""),
                                   db = "all", 
                                   result_dir = input_table[i, "ResultFileNames"],
                                   file_id = input_table[i, "File_id"], 
                                   input_dir, 
                                   no_of_candidates = 30,
                                   ppmx = 15)
    
    
    # Extract MS2 peak lists
    spec_pr2 <- ms2_peaks(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "."), sep =""), 
                          proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "."), sep =""),
                          input_dir,
                          result_dir = input_table[i, "ResultFileNames"],
                         file_id = input_table[i, "File_id"])
    
    # camera results for isotopes
    cam_res <- cam_func(input_dir,
                        f = input_table[i, "mzml_files"], 
                        ms2features = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/MS2DATA.csv", sep = ""), "."), sep = ""))
    
    # Extract MS1 peaks or isotopic peaks
    ms1p <- ms1_peaks(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS2DATA.csv', sep = ""), "."), sep =""), 
                      y = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/CAMERAResults.csv', sep = ""), "."), sep =""), 
                      input_table[i, "ResultFileNames"], 
                      input_dir, 
                      QCfile = TRUE)
    
    #prepare sirius parameter files
    sirius_param_files <- sirius_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA.csv', sep = ""), "."), sep =""), 
                                       result_dir = input_table[i, 'ResultFileNames'], 
                                       input_dir,
                                       SL = FALSE)
    
    
    # Run sirius
    run_sirius(files = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA_SiriusP.tsv', sep = ""), "."), sep =""), 
               ppm_max = 5, 
               ppm_max_ms2 = 15, 
               QC = TRUE, 
               SL = FALSE, 
               SL_path = NA,
               candidates = 30,
              profile = "qtof")
}

end.time <- Sys.time()

time.taken <- end.time - start.time
print(time.taken)
