#! /usr/bin/Rscript

#' @title Merge CAMERA results for multiple result csv files
#'
#' @description
#'
#' Exactly as the title says, perform this after the function
#' cam_fundMode


#' @param path where QC files are stored; store in a folder named QC
#'

#' @return
#' 
#' /QC/Combined_Camera_neg.csv
#' /QC/Combined_Camera_pos.csv
#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' merge_qc("/usr/project/QC/")


# ---------- Preparations ----------
# Load libraries
library("Spectra")
library("stringr")
library("dplyr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

path <- as.character(args[1])


# ---------- merge_qc ----------

merge_qc<- function(path){
    # combine all QC which are in positive mode
    df_pos <- list.files(path, pattern = "posCAMERA_Results_", full.names = TRUE) %>% 
        lapply(read_csv) %>% 
        bind_rows
    # remove any duplicated rows
    df_pos <- as.data.frame(df_pos[!duplicated(df_pos), ])

    #extract isotope column numbers, the numbers represent the group of isotope
    nm_p <- regmatches(df_pos[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",df_pos[, "isotopes"]))

    # for all the numbers, extract only first number, since it is the group number, 
    # second number can be charge
    for (i in 1:length(nm_p)){
        y <- as.numeric(unlist(nm_p[i]))
        df_pos[i,'istops'] = y[1]
    }

    # write csv for the combined_camera_pos results
    write.csv(df_pos, paste(path, "/Combined_Camera_pos.csv", sep = ""))
    
    # combine all QC which are in negative mode
    df_neg <- list.files(path, pattern = "negCAMERA_Results_", full.names = TRUE) %>% 
        lapply(read_csv) %>% 
        bind_rows
    # remove any duplicated rows based on mz
    df_neg <- as.data.frame(df_neg[!duplicated(df_neg), ])

    #extract isotope column numbers, the numbers represent the group of isotope
    nm_n <- regmatches(df_neg[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",df_neg[, "isotopes"]))

    # for all the numbers, extract only first number, since it is the group number, 
    # second number can be charge
    for (i in 1:length(nm_n)){
        y <- as.numeric(unlist(nm_n[i]))
        df_neg[i,'istops'] = y[1]
    }
    # write csv for the combined_camera_neg results
    write.csv(df_neg, paste(path, "/Combined_Camera_neg.csv", sep = ""))
}
merge_qc(path)
