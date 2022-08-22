start.time <- Sys.time()

# ---------- Script ----------
# input directory
input_dir <- getwd()
#input_dir <- "/Users/mahnoorzulfiqar/OneDriveUNI/MAW-data/MTBLS709_Data_32"
#input_dir

# load the functions file
source(paste(input_dir, "/Workflow_R_Functions.r", sep = ""))
#source("/Users/mahnoorzulfiqar/OneDriveUNI/MAW/Workflow_R_Functions.r")

input_table <- data.frame(ms2_rfilename(input_dir))
#input_table

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
    
}

end.time <- Sys.time()

time.taken <- end.time - start.time
print(time.taken)
