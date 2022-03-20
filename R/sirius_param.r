#! /usr/bin/Rscript

#' @title Write SIRIUS parameter files
#'
#' @description
#'
#' This function generates the parameter files for SIRIUS and saves 
#' the paths and names of the input and outputs in a csv along with
#' whether a isotope is present or not


#' @param x is the result csv file from the function - ms1peaks
#'
#' @param result_dir result directory for the current MS2 mzml file
#'
#' @param input_dir where all the input files are and the QC/ folder is
#'
#' @param SL either TRUE or FALSE (depends is the user has used 
#' inhouse library or here called a suspect list)
#'


#' @return
#' 
#' parameter .ms files for SIRIUS in a folder with the same name 
#' /insilico/SIRIUS/param1.ms
#' 
#' saves a csv file of the list of all the .ms parameter files and the
#' paths for the json result folders for Sirius results, another column
#' in this csv contains whether isotopic peaks were annotated or not
#'
#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' sirius_param(x = /usr/project/file1/insilico/MS1DATA.csv',
#'                      result_dir = "/usr/project/file1", 
#'                      input_dir = "/usr/project/" 
#'                      SL = TRUE)


# ---------- Preparations ----------
# Load libraries
library("stringr")
library("dplyr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

x <- as.character(args[1])
result_dir <- as.character(args[2])
input_dir <- as.character(args[3])
SL <- as.character(args[4])

# ---------- sirius_param ----------

# input x is result from either ms1_peaks
# SL is if a suspect list is present

sirius_param <- function(x, result_dir, input_dir, SL = TRUE){
    
    dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/SIRIUS", sep =""), "./"), sep = "")
    if (!file.exists(dir_name)){
        dir.create(dir_name, recursive = TRUE) ##create folder
    }
    isotopes <- c() #NA or isotope group number
    sirius_param_file <- c() #input for SIRIUS
    outputNames <- c() #output for SIRIUS with all db
    outputNamesSL <- c() #output for SIRIUS with suspect list as db

    parameter_file <- c()
    par <- 0
    
    x <- read.csv(x)
    
    a <-0 # counting
    y <- 0 # counting
    z <- 0 # counting
    for (i in 1:nrow(x)){
        
        par <- par+1
        para <- as.character(par) # for numbering
        
        #no MS1 PEAKS and no ISOTOPES
        
        if (x[i, "ms1Peaks"] == 'no ms1 peaks in QC'){

            #INPUT FILE NAME
            fileR <- paste(dir_name, "/" ,para, "_NA_iso_NA_MS1p_", x[i, "premz"], "_SIRIUS_param.ms", sep = "")
            
            sirius_param_file <- c(sirius_param_file, fileR)
            #ISOTOPE Information
            isotopes <- c(isotopes, NA)
            #OUTPUT
            fileSR <- paste(str_remove(fileR, ".ms"),'.json', sep = '')
            fileSRS <- paste(str_remove(fileR, ".ms"), 'SList.json', sep = '')
            outputNames <- c(outputNames, fileSR)
            outputNamesSL <- c(outputNamesSL, fileSRS)
            file.create(fileR, recursive = TRUE)
            file.conn <- file(fileR)
            open(file.conn, open = "at")
            
            #compound
            writeLines(paste(">compound", x[i,"id_X"], sep=" "),con=file.conn)
            #parentmass
            writeLines(paste(">parentmass", x[i,"premz"], sep=" "),con=file.conn)
            ##charge
            if (x[i,"pol"] == "pos"){
                writeLines(paste(">charge", "+1" ,sep=" "),con=file.conn)
            }
            else{
                writeLines(paste(">charge", "-1" ,sep=" "),con=file.conn)
            }
            #rt
            writeLines(paste(">rt", paste(x[i,"rtmed"], "s", sep =''), sep=" "),con=file.conn)
            
            #ms1
            writeLines(">ms1",con=file.conn)
            writeLines(paste(x[i,"premz"], x[i,"int"] ,sep=" "),con=file.conn)
            
            #ms2
            writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            
            ms2pk <- paste(input_dir, str_remove(x[i,"ms2Peaks"], "./"), sep ="")
            
            peak<- read.table(ms2pk)
            for (k in 1:length(peak[,1])){
                writeLines(paste(as.character(peak[k,1]),as.character(peak[k,2]), sep =" "), con=file.conn) 
            }
            close(file.conn)
            parameter_file <- c(parameter_file,file.conn)
        }
        
        # MS1 PEAKS and no ISOTOPES
        
        else if (grepl("_no_isotopes.txt", x[i, "ms1Peaks"], fixed=TRUE)){

            #INPUT FILE NAME
            fileR <- paste(dir_name, "/", para, "_NA_iso_MS1p_", x[i, "premz"], "_SIRIUS_param.ms", sep = "")
            sirius_param_file <- c(sirius_param_file, fileR)
            #ISOTOPE Information
            isotopes <- c(isotopes, NA)
            #OUTPUT
            fileSR <- paste(str_remove(fileR, ".ms"),'.json', sep = '')
            fileSRS <- paste(str_remove(fileR, ".ms"), 'SList.json', sep = '')
            outputNames <- c(outputNames, fileSR)
            outputNamesSL <- c(outputNamesSL, fileSRS)
            file.create(fileR, recursive = TRUE)
            file.conn <- file(fileR)
            open(file.conn, open = "at")
            
            #compound
            writeLines(paste(">compound", x[i,"id_X"], sep=" "),con=file.conn)
            #parentmass
            writeLines(paste(">parentmass", x[i,"premz"], sep=" "),con=file.conn)
            ##charge
            if (x[i,"pol"] == "pos"){
                writeLines(paste(">charge", "+1" ,sep=" "),con=file.conn)
            }
            else{
                writeLines(paste(">charge", "-1" ,sep=" "),con=file.conn)
            }
            #rt
            writeLines(paste(">rt", paste(x[i,"rtmed"], "s", sep =''), sep=" "),con=file.conn)
            
            #ms1
            writeLines(">ms1",con=file.conn)
            
            ms1pk <- paste(input_dir, str_remove(x[i,"ms1Peaks"], "./"), sep ="")
            peakms1<- read.table(ms1pk)
            
            for (l in 1:length(peakms1[,1])){
                writeLines(paste(as.character(peakms1[l,1]),as.character(peakms1[l,2]), sep =" "), con=file.conn) 
            }
            
            #ms2
            writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            
            ms2pk <- paste(input_dir, str_remove(x[i,"ms2Peaks"], "./"), sep ="")
            
            peakms2<- read.table(ms2pk)

            for (k in 1:length(peakms2[,1])){
                writeLines(paste(as.character(peakms2[k,1]),as.character(peakms2[k,2]), sep =" "), con=file.conn) 
            }
            
            close(file.conn)
            parameter_file <- c(parameter_file,file.conn)
        }
        
        # MS1 PEAKS and ISOTOPES
        
        else if (grepl("_isotopeNum_", x[i, "ms1Peaks"], fixed=TRUE)){

            #INPUT FILE NAME
            fileR <- paste(dir_name, "/", para, "_iso_MS1p_", as.character(x[i, "premz"]), "_SIRIUS_param.ms", sep = "")
            sirius_param_file <- c(sirius_param_file, fileR)
            #ISOTOPE Information
            isotopes <- c(isotopes, "present")
            #OUTPUT
            fileSR <- paste(str_remove(fileR, ".ms"),'.json', sep = '')
            fileSRS <- paste(str_remove(fileR, ".ms"), 'SList.json', sep = '')
            outputNames <- c(outputNames, fileSR)
            outputNamesSL <- c(outputNamesSL, fileSRS)
            file.create(fileR, recursive = TRUE)
            file.conn <- file(fileR)
            open(file.conn, open = "at")
            
             #compound
            writeLines(paste(">compound", x[i,"id_X"], sep=" "),con=file.conn)
            #parentmass
            writeLines(paste(">parentmass", x[i,"premz"], sep=" "),con=file.conn)
            ##charge
            if (x[i,"pol"] == "pos"){
                writeLines(paste(">charge", "+1" ,sep=" "),con=file.conn)
            }
            else{
                writeLines(paste(">charge", "-1" ,sep=" "),con=file.conn)
            }
            #rt
            writeLines(paste(">rt", paste(x[i,"rtmed"], "s", sep =''), sep=" "),con=file.conn)
            
            
            #ms1
            writeLines(">ms1",con=file.conn)
            
            ms1pk <- paste(input_dir, str_remove(x[i,"ms1Peaks"], "./"), sep ="")
            peakms1<- read.table(ms1pk)
            
            for (l in 1:length(peakms1[,1])){
                writeLines(paste(as.character(peakms1[l,1]),as.character(peakms1[l,2]), sep =" "), con=file.conn) 
            }
            
            #ms2
            writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            
            ms2pk <- paste(input_dir, str_remove(x[i,"ms2Peaks"], "./"), sep ="")
            
            peakms2<- read.table(ms2pk)

            for (k in 1:length(peakms2[,1])){
                writeLines(paste(as.character(peakms2[k,1]),as.character(peakms2[k,2]), sep =" "), con=file.conn) 
            }
            
            close(file.conn)
            parameter_file <- c(parameter_file,file.conn)
            
        }
    }
    if (SL){
        
        write.csv(data.frame(cbind(sirius_param_file, outputNames, outputNamesSL, isotopes)), paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA_SiriusPandSL.csv', sep = ""), "./"), sep =""))
        return(data.frame(cbind(sirius_param_file, outputNames, outputNamesSL, isotopes)))
        
    }
    else{
        
        write.csv(data.frame(cbind(sirius_param_file, outputNames, isotopes)), paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA_SiriusP.csv', sep = ""), "./"), sep =""))
        return(data.frame(cbind(sirius_param_file, outputNames, isotopes)))
        
    }
    
}


sirius_param(x, result_dir, input_dir, SL)

