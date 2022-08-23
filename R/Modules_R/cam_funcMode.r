#! /usr/bin/Rscript

#' @title CAMERA function with unknown mode
#'
#' @description
#'
#' This function is used to extract the modes of the spectra in the
#' QC files and separate them into pos and neg mode QC mzml files which
#' are then read by CAMERA and used for isotope annotations. The results 
#' from annotation are stored in csv file. When more files, use merge_qc 
#' after this function to merge the results of all QC in one csv to be 
#' used with all MS2 spectra files

#'
#' @param path where QC files are stored; store in a folder named QC
#'
#' @param pattern define this so that the function only catches the file 
#' with a certain pattern e.g: “common”, if no pattern, by default 
#' function takes “.mzML” as the pattern which is all the .mzML files
#'

#' @return
#' If single mode files:
#' CAMERA results in csv
#' If double modes:
#' Modes separated into different 
#' CAMERA results in csv

#' 
#' 

#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' cam_funcMode("/usr/project/QC/", ".mzmL")
#' 

# ---------- Preparations ----------
# Load libraries
library("Spectra")
library("stringr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)

path <- as.character(args[1])
pattern <- as.character(args[2])

# ---------- ms2_peaks ----------

# for QC files, which have both polarities in one file or the polarities are unknown
# inputs:
#path e.g: “/users/name/project/QC”
#pattern (define this so that the function only catches the file 
#with a certain pattern e.g: “common”, if no pattern, by default 
#function takes “.mzML” as the pattern which is all the .mzML files, 
#so make sure only the QC files are in this directory and not the MS2.mzML files)

cam_funcMode <- function(path, pattern = ".mzML"){
    library("CAMERA")
    # List all files present in QC folder 
    files_QC_N <- list.files(path, pattern = pattern ,full.names=TRUE)
    
    for (i in 1:length(files_QC_N)){
    
        # read each file using Spectra
        sps_all <- Spectra(files_QC_N[i], backend = MsBackendMzR())
        
        
        if (length(unique(sps_all$polarity)) == 1){
            if(unique(sps_all$polarity) == 1){
                #Read the same file with MS1 information; note CAMERA reads xcmsSet object
                xs <- xcmsSet(file = as.character(files_QC_N[i]),
                              profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                              mslevel= 1, progressCallback=NULL, polarity="positive",
                              scanrange = NULL, BPPARAM = bpparam(),
                              stopOnError = TRUE)
                # Create an xsAnnotate object 
                an <- xsAnnotate(xs) 
                # Group based on RT 
                anF <- groupFWHM(an, perfwhm = 0.6)
                # Annotate isotopes 
                anI <- findIsotopes(anF, mzabs = 0.01) 
                # Verify grouping 
                anIC <- groupCorr(anI, cor_eic_th = 0.75)
                #Annotate adducts 
                anFA <- findAdducts(anIC, polarity="positive") 
                #get a feature list
                peaklist <- getPeaklist(anFA)
                # add file_origin information
                peaklist$file_origin <- as.character(files_QC_N[i])
                file_name <- paste(path,"/posCAMERA_Results_", i, ".csv", sep = "")
                # write individual QC files which are in pos mode (later code will combine them)
                write.csv(peaklist, file = file_name)
            }else if (unique(sps_all$polarity) == 0){
                #Read the same file with MS1 information; note CAMERA reads xcmsSet object
                xs <- xcmsSet(file = as.character(files_QC_N[i]),
                              profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                              mslevel= 1, progressCallback=NULL, polarity="negative",
                              scanrange = NULL, BPPARAM = bpparam(),
                              stopOnError = TRUE)
                # Create an xsAnnotate object 
                an <- xsAnnotate(xs) 
                # Group based on RT 
                anF <- groupFWHM(an, perfwhm = 0.6)
                # Annotate isotopes 
                anI <- findIsotopes(anF, mzabs = 0.01) 
                # Verify grouping 
                anIC <- groupCorr(anI, cor_eic_th = 0.75)
                #Annotate adducts 
                anFA <- findAdducts(anIC, polarity="negative") 
                #get a feature list
                peaklist <- getPeaklist(anFA)
                # add file_origin information
                peaklist$file_origin <- as.character(files_QC_N[i])
                file_name <- paste(path, "/negCAMERA_Results_", i, ".csv", sep = "")
                # write individual QC files which are in pos mode (later code will combine them)
                write.csv(peaklist, file = file_name)
            }
        }
    else{
        pos <- sps_all[sps_all$polarity == 1]
        neg <- sps_all[sps_all$polarity == 0]
        file_p <- paste(path, "/QC_280k_pos", i, ".mzML", sep = "")
        file_n <- paste(path, "/QC_280k_neg", i, ".mzML", sep = "")
        
        # create new mzML QC files for pos and neg modes from each common QC file
        export(pos, backend = MsBackendMzR(), file = file_p)
        export(neg, backend = MsBackendMzR(), file = file_n)
        #Read the same file with MS1 information; note CAMERA reads xcmsSet object
        xs <- xcmsSet(file = as.character(file_p),
                        profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                        mslevel= 1, progressCallback=NULL, polarity="positive",
                        scanrange = NULL, BPPARAM = bpparam(),
                        stopOnError = TRUE)
        # Create an xsAnnotate object 
        an <- xsAnnotate(xs) 
        # Group based on RT 
        anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes 
        anI <- findIsotopes(anF, mzabs = 0.01) 
        # Verify grouping 
        anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts 
        anFA <- findAdducts(anIC, polarity="positive") 
        #get a feature list
        peaklist <- getPeaklist(anFA)
        # add file_origin information
        peaklist$file_origin <- as.character(files_QC_N[i])
        file_name <- paste(path, "/posCAMERA_Results_", i, ".csv", sep = "")
        # write individual QC files which are in pos mode (later code will combine them)
        write.csv(peaklist, file = file_name)
        
        #Read the same file with MS1 information; note CAMERA reads xcmsSet object
        xs <- xcmsSet(file = as.character(file_n),
                        profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                        mslevel= 1, progressCallback=NULL, polarity="negative",
                        scanrange = NULL, BPPARAM = bpparam(),
                        stopOnError = TRUE)
        # Create an xsAnnotate object 
        an <- xsAnnotate(xs) 
        # Group based on RT 
        anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes 
        anI <- findIsotopes(anF, mzabs = 0.01) 
        # Verify grouping 
        anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts 
        anFA <- findAdducts(anIC, polarity="negative") 
        #get a feature list
        peaklist <- getPeaklist(anFA)
        # add file_origin information
        peaklist$file_origin <- as.character(files_QC_N[i])
        file_name <- paste(path, "/negCAMERA_Results_", i, ".csv", sep = "")
        # write individual QC files which are in pos mode (later code will combine them)
        write.csv(peaklist, file = file_name)
    }
    
    }
    
    detach("package:CAMERA", unload=TRUE)
}

cam_funcMode(path, pattern)
