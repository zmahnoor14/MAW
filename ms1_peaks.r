#! /usr/bin/Rscript

#' @title Extract isotopic peaks and save as txt files
#'
#' @description
#'
#' This function takes the isoptpic peak annotations performed by CAMERA
#' and generates the MS1 or isotopic peak lists for each precursor m/z


#' @param x is the result csv file from the function - ms2peaks
#'
#' @param y is the combined camera results in csv format
#'
#' @param result_dir result directory for the current MS2 mzml file
#'
#' @param input_dir where all the input files are and the QC/ folder is
#'
#' @param QCfile either TRUE or FALSE (depends is the user has used QC or CAMERA)
#'


#' @return
#' 
#' a csv file containing all features and isotopic peak paths
#' isotopic peaks saved txt file for each precursor m/z
#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' ms1_peaks(x = /usr/project/file1/insilico/MS2DATA.csv',
#'                       y = "/usr/project/QC/combinedCam.csv", 
#'                      result_dir = "/usr/project/file1", 
#'                      input_dir = "/usr/project/" 
#'                      QC = TRUE)


# ---------- Preparations ----------
# Load libraries
library("stringr")
library("dplyr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

x <- as.character(args[1])
y <- as.character(args[2])
result_dir <- as.character(args[3])
input_dir <- as.character(args[4])
QCfile <- as.logical(args[5])

# ---------- ms1_peaks ----------

# Extract isotopic peaks for each pre_mz
# The input is x = first_list (from ms2peaks function) and y = camera results 

ms1_peaks <- function(x, y, result_dir, input_dir, QCfile = TRUE){
    # store the ms1_peak list path here
    ms1Peaks <- c()
    
    if (QCfile){
        
        dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/peakfiles_ms1", sep =""), "./"), sep = "")
        # create a new directory to store all the peak list txt files
        if (!file.exists(dir_name)){
            dir.create(dir_name, recursive = TRUE)
        }
        
        # read the CAMERA results
        y = read.csv(y)
        x = read.csv(x)
        
        # for all indices in the ms2 features table
        for (i in 1:nrow(x)){
            #store the indices of CAMERA that have same mz and rt as the ms2 features table
            store_c <- c()
            # for all indices in CAMERA Results
            for (j in 1:nrow(y)){
                # if mz and rt from ms2 features are within the range of CAMERA mz and rt
                if (x[i, 'premz'] <= y[j, "mzmax"]   && y[j, "mzmin"] <= x[i, 'premz'] && x[i, 'rtmed'] <= y[j, "rtmax"] && y[j, "rtmin"] <= x[i, 'rtmed']){
                    store_c <- c(store_c, j)
                }
            }
            # indices with same pre m/z and same rt
            df_y <- y[store_c, ]
            df_y <- as.data.frame(df_y)

            #if there was only one index
            if (nrow(df_y) == 1){
                
                # if there was no isotope annotation for that one index
                if (is.na(df_y[1, "istops"])){
                
                    mz <- df_y[1, "mz"] # save mz
                    int <- df_y[1, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(dir_name, "/ms1_peaks_", x[i, 'premz'], "_no_isotopes.txt", sep = "") # save name of the peaklist
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE) # save peak list
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1) # add the path of the peak list to a list
                }
                
                # if there was an isotope annotation
                else{
                
                    df_x <- y[which(y[, "file_origin"] ==df_y[1, "file_origin"]), ] # extract camera results from one file origin
                    df_x <- df_x[which(df_x[, 'istops'] == df_y[1, 'istops']), ] # extract only certain isotope annotation group
                    mz <- df_x[, "mz"] # save mz
                    int <- df_x[, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(dir_name, "/ms1_peaksISOTOPE_", x[i, 'premz'], "_isotopeNum_", df_x[1, "istops"], ".txt", sep = "")
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE)
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1)
                }
            }
            # if there are more indices for df_y
            else if(nrow(df_y) > 1){
                # if all enteries have no isotope annotation
                if(all(is.na(df_y[, 'istops']))){
                
                    df_z <- df_y[which(df_y[,"into"] == max(df_y[,"into"])), ] # extract the ms1 peak with highest intensity
                    mz <- df_z[1, "mz"] # save mz
                    int <- df_z[1, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(dir_name, "/ms1_peaks_", x[i, 'premz'], "_no_isotopes.txt", sep = "") # save name of the peaklist
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE) # save peak list
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1) # add the path of the peak list to a list
                }
                # if not all isotope annotations are NA
                else if (!(all(is.na(df_y[, 'istops'])))){
                
                    df_y <- df_y[!is.na(df_y$'istops'),] # Remove the NA isotope annotations
                    df_z <- df_y[which(df_y[,"into"] == max(df_y[,"into"])), ] # Select the MS1 peak with highest intensity
                    df_z1 <- y[which(y[, "file_origin"] == df_z[1, "file_origin"]), ]  # extract camera results from one file origin
                    df_z1 <- df_z1[which(df_z1[, 'istops'] == df_z[1, 'istops']), ] # extract only certain isotope annotation group
                    mz <- df_z1[, "mz"] # save mz
                    int <- df_z1[, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(dir_name, "/ms1_peaksISOTOPE_", x[i, 'premz'], "_isotopeNum_", df_z1[1, 'istops'],".txt", sep = "") # save name of the peaklist
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE) # save peak list
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1) # add the path of the peak list to a list
                }
            }
            else if (nrow(df_y)==0){
                ms1Peaks <- c(ms1Peaks, 'no ms1 peaks in QC')
            }
        }
        second_list <- data.frame(cbind(x, ms1Peaks))
        write.csv(second_list, file = paste(result_dir,'/insilico/MS1DATA.csv', sep = ""))
        return(second_list)
    }
    else{
        ms1Peaks <- c(ms1Peaks, 'no ms1 peaks in QC')
        second_list <- data.frame(cbind(x, ms1Peaks))
        write.csv(second_list, file = paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA.csv', sep = ""), "./"), sep =""))
        return(second_list)
    }
    
}

# Usage
ms1_peaks(x, y, result_dir, input_dir, QCfile)

