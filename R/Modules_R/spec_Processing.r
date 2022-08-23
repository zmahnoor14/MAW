#! /usr/bin/Rscript

#' @title Read mzML files and extract precursor m/z(s)
#'
#' @description
#'
#' Spectra package reads files in .mzML format, filtering any empty
#' spectrum and ommitting NA values from precurosr m/z list. 
#' 


#' @param input_dir is full directory where all MZML input files
#'
#' @param db is either one of the spectral libraries which can be 
#'        gnps, hmdb, mbank or all


#' @return
#' 
#' This function returns the spectra object and a list of 
#' precursor m/z(s) present in the .mzML file
#' 
#' 
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' spec_Processing("/usr/project/MS2specfile1.mzML", "/usr/project/MS2specfile1")


# ---------- Preparations ----------
# Load libraries
library("Spectra")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

x <- as.character(args[1])
result_dir <- as.character(args[2])

# ---------- spec_Processing ----------

#' All spectra in mzML files preprocessing, return two outputs, pre-processed MS2 spectra and all precursor masses
# x is one mzML file
#' All spectra in mzML files preprocessing, return two outputs, pre-processed MS2 spectra and all precursor masses
# x is one mzML file
spec_Processing <- function(input_dir, x, result_dir){
    
    input_dir <- dirname(arg[1])
    
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

