#! /usr/bin/Rscript

#' @title Run Metfrag 
#'
#' @description
#'
#' This function runs Metfrag and generates result csv files
#' 

#' @param met_param is the txt file conatining list of paths to metfrag parameter files
#' 
#' @param input_dir is the directory with all input files and necessary files stored
#'

#' @return
#' 
#' a csv result files for all metfrag parameter files
#' 
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' met_param <- metfrag_param(met_param = "usr/project/file1/insilico/metparam_list.txt",
#'                                input_dir = "usr/project/")

# ---------- Preparations ----------
# Load libraries

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

met_param <- as.character(args[1])
input_dir <- as.character(args[2])


# ---------- run_metfrag ----------

run_metfrag <- function(met_param, input_dir){
    
    filesmet_param <- read.table(met_param)
    
    for (files in filesmet_param[[1]]){
        system(paste("java -jar",  paste(input_dir, "MetFrag2.4.5-CL.jar", sep = ''), files))
    }
}
run_metfrag(met_param, input_dir)
