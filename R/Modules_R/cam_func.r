#! /usr/bin/Rscript

#' @title Normal CAMERA
#'
#' @description
#'
#' CAMERA is used to extract isotopic peaks from the qc files
#' depending on the mode of the single QC.mzml file


#' @param path where QC files are stored; store in a folder named QC
#'
#' @param f file name
#'
#' @param mode either type pos or neg depending on mode of the qc.mzml
#'
#' @param input_dir where all the input files are and the QC/ folder is
#'

#' @return
#' 
#' /QC/negCAMERAResults_qc1.csv
#' /QC/posCAMERAResults_qc1.csv
#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' cam_func("/usr/project/QC/", "QC1pos.mzML", "pos", "/usr/project/")


# ---------- Preparations ----------
# Load libraries
library("stringr")
library("dplyr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

path <- as.character(args[1])
f <- as.character(args[2])
mode <- as.character(args[3])
input_dir <- as.character(args[4])

# ---------- came_func ----------

cam_func <- function(path, f, mode = "pos", input_dir){
    library("CAMERA")
    fl <- paste(path, f, sep ="")
    if(mode == "pos"){
        xs <- xcmsSet(file = fl,profmethod = "bin", 
              profparam = list(), lockMassFreq=FALSE,
              mslevel= 1, progressCallback=NULL, polarity="positive",
              scanrange = NULL, BPPARAM = bpparam(),stopOnError = TRUE)
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
        peaklist <- getPeaklist(anFA)
        peaklist$file_origin <- fl

        #extract isotope column numbers, the numbers represent the group of isotope
        nm_po <- regmatches(peaklist[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",peaklist[, "isotopes"]))
        # for all the numbers in v, extract only first number, since it is the group number, 
        # second number can be charge
        for (i in 1:length(nm_po)){
            y <- as.numeric(unlist(nm_po[i]))
            peaklist[i,'istops'] = y[1]
        }
        name <- str_remove(f, ".mzML")
        write.csv(peaklist, file = paste(input_dir, "QC/posCAMERAResults_", name,".csv", sep = ""))
    }
    else if(mode == "neg"){
        xs <- xcmsSet(file = fl,profmethod = "bin", 
              profparam = list(), lockMassFreq=FALSE,
              mslevel= 1, progressCallback=NULL, polarity="negative",
              scanrange = NULL, BPPARAM = bpparam(),stopOnError = TRUE)
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
        peaklist <- getPeaklist(anFA)
        peaklist$file_origin <- fl

        #extract isotope column numbers, the numbers represent the group of isotope
        nm_ne <- regmatches(peaklist[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",peaklist[, "isotopes"]))
        # for all the numbers in v, extract only first number, since it is the group number, 
        # second number can be charge
        for (i in 1:length(nm_ne)){
            y <- as.numeric(unlist(nm_ne[i]))
            peaklist[i,'istops'] = y[1]
        }
        name <- str_remove(f, ".mzML")
        write.csv(peaklist, file = paste(input_dir, "QC/negCAMERAResults_", name,".csv", sep = ""))
    }
    detach("package:CAMERA", unload=TRUE)
    
}
cam_func(path, f, mode, input_dir)
