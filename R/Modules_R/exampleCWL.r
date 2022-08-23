#! /usr/bin/Rscript

library(Spectra)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
#print(args)

input_dir <- as.character(args[1])

## Specifying a function for creating result directories for each input mzml
# input for the function:
# input directory
ms2_rfilename<- function(input_dir){
    
    input_dir <- dirname(input_dir)
    
    if (dir.exists(input_dir)){
        #list_ms2_files <- intersect(list.files(input_dir, pattern = "_PRM_"), list.files(input_dir, pattern = ".mzML"))
        list_ms2_files <- list.files(input_dir, pattern = ".mzML")
        mzml_file <- paste(input_dir, "/", list_ms2_files, sep = "")

        #store the result file names to return to this function as output
        mzml_files <- c()
        ResultFileNames <- c()
        File_id <- c()
        nx <- 0
        # x is mzML files
        for (i in 1:length(mzml_file)){
            nx <- nx+1
            # remove .mzML to extract just the names
            mzml_filex <- str_replace(mzml_file[i], input_dir, ".")
            name_mzmls <- str_remove(as.character(mzml_filex), ".mzML")
            name_mzmlsd <- str_remove(mzml_file[i], ".mzML")
            if (!file.exists(name_mzmlsd)){
                dir.create(name_mzmlsd) ##create folder
            }
            ResultFileNames<- c(ResultFileNames, name_mzmls)
            mzml_files <- c(mzml_files, mzml_filex)
            File_id <- c(File_id, paste("file_", nx, sep = ""))
        }
        input_table <- cbind(mzml_files, ResultFileNames, File_id)

        write.csv(input_table, paste(input_dir, "/input_table.csv", sep = ""))
        return(data.frame(input_table))
    }
    else{
        stop("Your input_dir is incorrect. Please provide the directory where all your input files are stored. An example would be: '/Users/my_name/input_dir/'. don't forget the '/' at the end : ) Good Luck")
    }
}

spec_Processing <- function(input_dir, x, result_dir){
    
    input_dir <- dirname(input_dir)
    
    x <- paste(input_dir, str_remove(x, "."), sep = "")
    
    result_dir <- paste(input_dir, str_remove(result_dir, "."), sep = "")
    
    if (file.exists(x) && substring(x, nchar(x)) == "L"){
        if (dir.exists(result_dir)){
            # read the spectra
            sps_all <- Spectra(x, backend = MsBackendMzR())
            #' Change backend to a MsBackendDataFrame: load data into memory
            #sps_all <- setBackend(sps_all, MsBackendDataFrame())
            #' Filter Empty Spectra
            sps_all <- filterEmptySpectra(sps_all)
            #' Extract Precursor m/z(s) in each file
            pre_mz <- unique(precursorMz(sps_all))
            #' Remove any NAs
            pre_mz <- na.omit(pre_mz)
            export(sps_all, backend = MsBackendMzR(), file = paste(result_dir, "/processedSpectra.mzML", sep = ""))
            write.table(pre_mz, file = paste(result_dir, "/premz_list.txt", sep = ""), sep = "/t", row.names = FALSE, col.names = FALSE)
            spsall_pmz <- list(sps_all, pre_mz)
            return(spsall_pmz)
        }
        else{
            stop("Seems like it is not the result directory of the input .mzML file which is provided as x. Please use the function ms2_rfilename to generate a result directory or create one yourself with the same name as the .mzML input file.")
        }
    }
    else{
        stop("Are you sure x is an mzML input file?")
    }
    
}

input_table <- data.frame(ms2_rfilename(input_dir))

spec_pr <- spec_Processing(input_dir,
                           input_table[1, "mzml_files"], 
                           input_table[1, "ResultFileNames"])
