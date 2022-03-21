#! /usr/bin/Rscript

#' @title Metfrag parameter files
#'
#' @description
#'
#' This function is used to generate the txt parameter file for MetFrag
#' 

#' @param x is the csv file resulting from sirius_postprocessing
#' 
#' @param result_dir  is the result directory for the input MS2 mzml file
#' 
#' @param input_dir is the directory with all input files and necessary files stored
#' 
#' @param adducts is the csv file that conatains information on adducts notation used in Metfrag
#' 
#' @param sl_mtfrag is the txt file of inchiKeys which is used by MetFrag
#'
#' @param SL either TRUE or FALSE (depends is the user has used 
#' inhouse library or here called a suspect list)
#' 
#' @param ppm_max the ppm range to be explored with the precursor masses found in the DBs
#'
#' @param ppm_max_ms2, the ppm range to be explored with the fragmenat peaks found in the DBs
#'

#' @return
#' 
#' a txt file containing lis of paths of metfrag parameter files
#' 
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' met_param <- metfrag_param(x = "usr/project/file1/insilico/MS1DATAsirius.csv",
#'                                result_dir = "usr/project/file1",
#'                                input_dir = "usr/project/",
#'                                adducts = "usr/project/MetFrag_AdductTypes.csv", 
#'                                sl_mtfrag = "usr/project/sl_metfrag.txt", 
#'                                SL = TRUE,
#'                                ppm_max = 5, 
#'                                ppm_max_ms2= 15)


# ---------- Preparations ----------
# Load libraries
library(stringr)
library(readr)
library(dplyr)

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

x <- as.character(args[1])
result_dir <- as.character(args[2])
input_dir <- as.character(args[3])
adducts <- as.character(args[4])
sl_mtfrag <- as.character(args[5])
SL <- as.logical(args[6])
ppm_max <- as.numeric(args[7])
ppm_max_ms2 <- as.numeric(args[8])


# ---------- metfrag_param ----------

#x is the dataframe result from sirius_postprocess
#result_dir 
#input_dir 
#adducts 
#sl_mtfrag 
#SL = TRUE if suspect list present

metfrag_param <- function(x, result_dir, input_dir, adducts, sl_mtfrag, SL = TRUE, ppm_max = 5, ppm_max_ms2= 15){
    
    x <- read.csv(x)
    
    dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/MetFrag", sep =""), "./"), sep = "")
    
    if (!file.exists(dir_name)){
        dir.create(dir_name, recursive = TRUE) ##create folder
    }
    
    db <- c("PubChem", "KEGG")
    
    parameter_file <- c()
    par <- 0
    metfrag_param_file <- c()
    
    AdductsMF <- data.frame(read.csv(adducts))
    
    for (j in 1:nrow(x)){
        if (!(is.na(x[j, 'Adducts']))){
            for (k in db){
                par <- par+1
                para <- as.character(par)
                fileR <- paste(dir_name, "/", para, "_id_", x[j, 'id_X'], "_mz_", x[j, 'premz'], "_rt_", x[j, 'rtmed'], "_db_", k, ".txt", sep = '')
                metfrag_param_file <- c(metfrag_param_file, fileR)

                file.create(fileR, recursive = TRUE)
                file.conn <- file(fileR)
                open(file.conn, open = "at")
                
                
                peakspath <- str_replace(x[j, "ms2Peaks"], "./", input_dir)
                
                
                #writeLines(paste("PeakListPath = ",as.character(peakspath),sep=""),con=file.conn)
                writeLines(paste("PeakListPath = ",peakspath, sep=""),con=file.conn)
                writeLines(paste("IonizedPrecursorMass = ", x[j, "premz"]),con = file.conn)
                
                # write code here
                
                
                PrecursorIonMode <- AdductsMF[which(AdductsMF[, "PrecursorIonType"] == gsub("[[:blank:]]", "", x[j, 'Adducts'])), "PrecursorIonMode"]
                IsPositiveIonMode <- AdductsMF[which(AdductsMF[, "PrecursorIonType"] == gsub("[[:blank:]]", "", x[j, 'Adducts'])), "IsPositiveIonMode"]
            
                
                writeLines(paste("PrecursorIonMode = ", PrecursorIonMode, sep = ''), con = file.conn)
                writeLines(paste("IsPositiveIonMode = ", IsPositiveIonMode, sep = ''), con = file.conn)
                
                writeLines(paste("MetFragDatabaseType = ", k),con = file.conn)
            
                writeLines(paste("DatabaseSearchRelativeMassDeviation = ", ppm_max, sep = ''),con=file.conn)
                writeLines("FragmentPeakMatchAbsoluteMassDeviation = 0.001",con=file.conn)
                writeLines(paste("FragmentPeakMatchRelativeMassDeviation = ", ppm_max_ms2, sep = ''),con=file.conn)
                
                if (SL){
                    writeLines(paste("ScoreSuspectLists = ", sl_mtfrag),con=file.conn)
            
                    writeLines("MetFragScoreTypes = FragmenterScore, SuspectListScore",con=file.conn)
                    writeLines("MetFragScoreWeights = 1.0, 1.0", con=file.conn)
                }
                else{
                    writeLines("MetFragScoreTypes = FragmenterScore", con=file.conn)
                    writeLines("MetFragScoreWeights = 1.0", con=file.conn)
                }
                
                
                writeLines("MetFragCandidateWriter = CSV",con=file.conn)
            
                writeLines(paste("SampleName = ", para, "_id_", x[j, 'id_X'], "_mz_", x[j, 'premz'], "_rt_", x[j, 'rtmed'], "_db_", k, sep = ''),con=file.conn)
                resultspath <- str_replace(result_dir, "./", input_dir)
                writeLines(paste("ResultsPath = ", resultspath, "/insilico/MetFrag/", sep = ''),con=file.conn)
            
                writeLines("MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter",con=file.conn)
                writeLines("MetFragPostProcessingCandidateFilter = InChIKeyFilter",con=file.conn)
                writeLines("MaximumTreeDepth = 2",con=file.conn)
                writeLines("NumberThreads = 1",con=file.conn)
                
                close(file.conn)
                parameter_file <- c(parameter_file,file.conn)
            
            }
        }
    }
    
    write.table(metfrag_param_file, file = paste(input_dir, str_remove(paste(result_dir, "/insilico/metparam_list.txt", sep =""), "./"), sep = ""), sep = "/t", row.names = FALSE, col.names = FALSE)
    return(metfrag_param_file)
}

# Usage:
# metfrag_param(x, result_dir, input_dir, adducts, sl_mtfrag, SL = TRUE)

metfrag_param(x, result_dir, input_dir, adducts, sl_mtfrag, SL, ppm_max, ppm_max_ms2)
