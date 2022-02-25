#! /usr/bin/Rscript

#' @title Result Directories 
#'
#' @description
#'
#' This function creates a dataframe of all the input files and their result 
#' directories (which are also created with this function). The function gives 
#' each input file a file id.
#' 


#' @param input_dir 


#' @return
#' 
#' A dataframe of mzml_files, ResultFileNames, File_id columns. 
#' Result files for each input .mzML is generated.
#' 
#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' ms2_rfilename("usr/project/")




# ---------- Preparations ----------
# Load libraries
library("stringr")



# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

input_dir <- as.character(args[1])


# ---------- ms2_rfilename ----------

## Specifying a function for creating result directories for each input mzml
# input for the function:
# input directory
ms2_rfilename<- function(input_dir){
    #list_ms2_files <- intersect(list.files(input_dir, pattern = "_PRM_"), list.files(input_dir, pattern = ".mzML"))
    list_ms2_files <- list.files(input_dir, pattern = ".mzML")
    mzml_file <- paste(input_dir, list_ms2_files, sep = "")
    
    #store the result file names to return to this function as output
    mzml_files <- c()
    ResultFileNames <- c()
    File_id <- c()
    nx <- 0
    # x is mzML files
    for (i in 1:length(mzml_file)){
        nx <- nx+1
        # remove .mzML to extract just the names
        mzml_filex <- str_replace(mzml_file[i], input_dir, "./")
        name_mzmls <- str_remove(as.character(mzml_filex), ".mzML")
        #name_mzml <- str_replace(name_mzmls, input_dir, "./")
        #' for each file a subdirectory is created to store all results in that, add working directory
        if (!file.exists(name_mzmls)){
            dir.create(name_mzmls) ##create folder
        }
        ResultFileNames<- c(ResultFileNames, name_mzmls)
        mzml_files <- c(mzml_files, mzml_filex)
        File_id <- c(File_id, paste("file_", nx, sep = ""))
    }
    input_table <- cbind(mzml_files, ResultFileNames, File_id)
    
    write.csv(input_table, paste(input_dir, "input_table.csv", sep = ""))
    return(data.frame(input_table))
}

ms2_rfilename(input_dir)
