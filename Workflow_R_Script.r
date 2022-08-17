#packageVersion("Spectra")
#packageVersion("MsBackendMgf")
#packageVersion("MsBackendHmdb")
#packageVersion("MsBackendMsp")
#packageVersion("MsCoreUtils")
#packageVersion("readr")
#packageVersion("dplyr")
#packageVersion("rvest")
#packageVersion("stringr")
#packageVersion("xml2")
#packageVersion("CAMERA")



# ---------- Script ----------
# input directory
input_dir <- getwd()
input_dir

# load the functions file
source(paste(input_dir, "/Workflow_R_Functions.r", sep = ""))

# Track Time 
#start_time <- Sys.time()

line <- system.time({
    # Run the first function; this creates a dataframe of your input files, their result directories 
    # and gives an id to each input file; stores the table in directory as a csv filr
    input_table <- data.frame(ms2_rfilename(input_dir))
    input_table
    for (i in 1:nrow(input_table)){

        #Preprocess and Read the mzMLfiles
        spec_pr <- spec_Processing(input_dir,
                                   input_table[i, "mzml_files"], 
                                   input_table[i, "ResultFileNames"])


        #perform dereplication with all dbs
        df_derep <- spec_dereplication(mzml_file = input_table[i, "mzml_files"],
                                       pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "."), sep =""), 
                                       proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "."), sep =""),
                                       db = "all", 
                                       result_dir = input_table[i, "ResultFileNames"],
                                       file_id = input_table[i, "File_id"], 
                                       input_dir, 
                                       ppmx = 15)

        # Extract MS2 peak lists
        spec_pr2 <- ms2_peaks(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "."), sep =""), 
                              proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "."), sep =""),
                              input_dir,
                              input_table[i, "ResultFileNames"],
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
<<<<<<< HEAD
                                           SL = FALSE)
=======
                                           SL = TRUE)
>>>>>>> refs/remotes/origin/main

        # Run sirius
        run_sirius(files = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA_SiriusPandSL.csv', sep = ""), "."), sep =""), 
                   ppm_max = 5, 
                   ppm_max_ms2 = 15, 
                   QC = TRUE, 
<<<<<<< HEAD
                   SL = FALSE, 
                   SL_path = NA,
=======
                   SL = TRUE, 
                   SL_path = paste(input_dir, 'SL_Frag/', sep = ""),
>>>>>>> refs/remotes/origin/main
                   candidates = 30)


        # Post process Sirius results and extract adducts for MetFrag
        adducts <- sirius_adduct(input_dir,
                                 x = input_table[i, "ResultFileNames"], 
<<<<<<< HEAD
                                 SL = FALSE)
=======
                                 SL = TRUE)
>>>>>>> refs/remotes/origin/main


        # prepare Metfrag parameter files
        met_param <- metfrag_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/MS1DATAsirius.csv", sep = ""), "."), sep =""), 
                                   result_dir = input_table[i, "ResultFileNames"],
                                   input_dir, 
<<<<<<< HEAD
                                   sl_mtfrag = NA, 
                                   SL = FALSE,
=======
                                   sl_mtfrag = paste(input_dir, "/SL_metfrag.txt", sep = ""), 
                                   SL = TRUE,
>>>>>>> refs/remotes/origin/main
                                   ppm_max = 5, 
                                   ppm_max_ms2= 15)


        # run metfrag
        run_metfrag(met_param = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/metparam_list.txt", sep = ""), "."), sep =""))


    }
})


write(paste("MAW took", line[1] "of seconds.", sep = "") ,file=paste(input_dir, "/summaryFile.txt", sep = ""),append=TRUE)

<<<<<<< HEAD
save.image(file = paste(input, "/StandardsRun1.Rdata", sep ="")
=======
save.image(file = paste(input, "/BryophytesRun1.Rdata", sep ="")
>>>>>>> refs/remotes/origin/main


