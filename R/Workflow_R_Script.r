library(parallel)
library(doParallel)
library(future)
library(iterators)
library(listenv)

options(future.globals.maxSize = 8 * 1024^3) # increase dataset size limit taken by future to 8GB

# detects number of cores
n.cores <- parallel::detectCores()


plan(list(
  tweak(multisession, workers = ((n.cores + 5) %/% 3) %/% 2),
  tweak(multisession, workers = ((n.cores + 5) %/% 3) %/% 2),
  tweak(multisession, workers = 3)
))


start.time <- Sys.time()

# ---------- Script ----------
# input directory
input_dir <- paste(getwd(), "/data", sep = "")
input_dir

# load the functions file
source(paste(getwd(), "/Workflow_R_Functions.r", sep = ""))

input_table <- data.frame(ms2_rfilename(input_dir))
input_table

input_table_idxs <- listenv()

for (i in 1:nrow(input_table)){

input_table_idxs[[i]] <- future({
    ##Preprocess and Read the mzMLfiles
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
                                            no_of_candidates = 50,
                                            ppmx = 15)
        Extract MS2 peak lists
        spec_pr2 <- ms2_peaks(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "."), sep =""),
                              proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "."), sep =""),
                              input_dir,
                              result_dir = input_table[i, "ResultFileNames"],
                             file_id = input_table[i, "File_id"])

        # camera results for isotopes
#         cam_res <- cam_func(input_dir,
#                             f = input_table[i, "mzml_files"], 
#                             ms2features = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/MS2DATA.csv", sep = ""), "."), sep = ""))

        # Extract MS1 peaks or isotopic peaks
        ms1p <- ms1_peaks(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS2DATA.csv', sep = ""), "."), sep =""),
                          y = NA, 
                          input_table[i, "ResultFileNames"],
                          input_dir,
                          QCfile = FALSE)


        #prepare sirius parameter files
        sirius_param_files <- sirius_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA.csv', sep = ""), "."), sep =""),
                                           result_dir = input_table[i, 'ResultFileNames'],
                                           input_dir,
                                           SL = FALSE)
        # Run sirius
        run_sirius(files = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA_SiriusP.tsv', sep = ""), "."), sep =""),
                   ppm_max = 5,
                   ppm_max_ms2 = 15,
                   QC = FALSE,
                   SL = FALSE,
                   SL_path = NA,
                   candidates = 30,
                  profile = "orbitrap",
                  db = "ALL")
   }) #end input_table_idxs future
}

input_table_idxs <- as.list(input_table_idxs)
v_input_table_idxs <- future::value(input_table_idxs)

end.time <- Sys.time()

time.taken <- end.time - start.time
print(time.taken)
