#!/Users/mahnoorzulfiqar/opt/anaconda3/envs/mawRpy/bin/Rscript

#' @title Result Directories
#'
#' @description
#'
#' It lists all files with .mzML extension.
#' These files are without the full path, so the input_dir is added to 
#' the path to these files. For each .mzML file, the extension is removed 
#' and stored as a name for the results directory for each input .mzML 
#' file For each input .mzML file an id is generated like file_1
#' A table is generated at the end


#' @param input_dir is full directory where all MZML input files
#'



#' @return
#' 
#' A table of mzml_files, ResultFileNames, File_id columns. 
#' Result files for each input .mzML is generated.
#' 

#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' ms2_rfilename(input_dir = "/usr/project/")


# ---------- Preparations ----------

# Load Libraries
library(stringr)


# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)


input_dir <- as.character(args[1])


## Specifying a function for creating result directories for each input mzml
# input for the function:
# input directory
ms2_rfilename<- function(input_dir){
    if (dir.exists(input_dir)){
        #list_ms2_files <- intersect(list.files(input_dir, pattern = "_PRM_"), list.files(input_dir, pattern = ".mzML"))
        list_ms2_files <- list.files(input_dir, pattern = ".mzML")
        mzml_file <- paste(input_dir, "/", list_ms2_files, sep = "")

        #store the result file names to return to this function as output
        mzml_files <- c()
        ResultFileNames <- c()
        File_id <- c()
        nx <- 0
        # x is mzML files
        for (i in 1:length(mzml_file)){
            nx <- nx+1
            # remove .mzML to extract just the names
            mzml_filex <- str_replace(mzml_file[i], input_dir, ".")
            name_mzmls <- str_remove(as.character(mzml_filex), ".mzML")
            name_mzmlsd <- str_remove(mzml_file[i], ".mzML")
            if (!file.exists(name_mzmlsd)){
                dir.create(name_mzmlsd) ##create folder
            }
            ResultFileNames<- c(ResultFileNames, name_mzmls)
            mzml_files <- c(mzml_files, mzml_filex)
            File_id <- c(File_id, paste("file_", nx, sep = ""))
        }
        input_table <- cbind(mzml_files, ResultFileNames, File_id)

        write.csv(input_table, paste(input_dir, "/input_table.csv", sep = ""))
        return(data.frame(input_table))
    }
    else{
<<<<<<< HEAD
        stop("Your input_dir is incorrect. Please provide the directory where all your input files are stored. An example would be: '/Users/my_name/input_dir'. Good Luck")
=======
        stop("Your input_dir is incorrect. Please provide the directory where all your input files are stored. An example would be: '/Users/my_name/input_dir/'. don't forget the '/' at the end : ) Good Luck")
>>>>>>> refs/remotes/origin/main
    }
}

ms2_rfilename(input_dir)
