#! /usr/bin/Rscript

#' @title CAMERA function with unknown mode
#'
#' @description
#'
#' for QC files, which have both polarities in one file or the polarities 
#' are unknown, use Spectra to divide them based on + and - polarities 
#' and CAMERA to annotate the isotopic peaks
#' 
#' @param path e.g: “/users/name/project/QC”
#' @param pattern (define this so that the function only catches the file 
#' with a certain pattern e.g: “common”, if no pattern, by default 
#' function takes “.mzML” as the pattern which is all the .mzML files, 
#' so make sure only the QC files are in this directory and not the MS2.mzML files)
#' 
#' 
#' @return
#' 
#' Modes separated into different .mzML files and CAMERA results in csv
#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' cam_funcMode("/usr/project/QC”, pattern = "OrbitrapMS2")

# ---------- Preparations ----------
# Load libraries
library("Spectra")
library("stringr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

path <- as.character(args[1])
pattern <- as.character(args[2])

# ---------- cam_funcMode ----------

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
