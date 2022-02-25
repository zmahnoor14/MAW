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
#' spec_Processing("usr/project/MS2specfile1.mzML")


# ---------- Preparations ----------
# Load libraries
library("Spectra")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

x <- as.character(args[1])


# ---------- spec_Processing ----------

#' All spectra in mzML files preprocessing, return two outputs, pre-processed MS2 spectra and all precursor masses
# x is one mzML file
spec_Processing <- function(x){
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
    spsall_pmz <- list(sps_all, pre_mz)
    return(spsall_pmz)
}
spec_Processing(x)
