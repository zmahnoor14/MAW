#! /usr/bin/Rscript

#' @title Download Spectral DBs as Spectra
#'
#' @description
#'
#' A summary file in a text format is saved in the input directory 
#' which stores the timing when each database is stored, e.g: GNPS is 
#' updated when there is a submission. This information is stored for 
#' reproducibility checks and versions of the databases.
#' If all is selected, then GNPS, HMDB, and MassBank all 
#' are downloaded with their URLs. If a specific database is selected only
#' that database is downloaded. Each database is stored as in the format 
#' it is downloadable from the database webpages.


#' @param input_dir is full directory where all MZML input files
#'
#' @param db is either one of the spectral libraries which can be 
#'        gnps, hmdb, mbank or all


#' @return
#' 
#' Spectral DB saved with the following name in the input_dir: 
#' gnps.rda, hmdb.rda, mbank.rda
#' summary file saved as summaryFile.txt which contains the 
#' timings and versions if available of saved databases
#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' download_specDB(input_dir = "/usr/project/", db = "all")


# ---------- Preparations ----------
# Load libraries
library("Spectra")
library("MsBackendMgf")
library("MsBackendMsp")
library("MsBackendHmdb")
library("rvest")
library("stringr")
library("xml2")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

input_dir <- as.character(args[1])
db <- as.character(args[2])
error <- as.logical(args[3])


# ---------- download_specDB ----------

download_specDB <- function(input_dir, db = "all", error = TRUE){
    # Track Time 
    start_time <- Sys.time()

    # only input available as of now
    databases <- 'gnps, hmdb, mbank, all'
    
    # creat a summary file, open and store timings of download and version if possible
    summaryFile <- paste(input_dir, "summaryFile.txt", sep = "")
    file.create(summaryFile, recursive = TRUE)
    file.conn <- file(summaryFile)
    open(file.conn, open = "at")
            
    # gnps
    if (db == "all" || db =="gnps"){
        
        print("GNPS WORKS")
        
        # Download file
        system(paste("wget -P", 
                     input_dir,
                     "https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.mgf",
                     sep =  " "))
        
        # load the spectra into MsBackendMgf
        gnpsdb <- Spectra(paste(input_dir, "ALL_GNPS.mgf", sep = ''), source = MsBackendMgf())
        save(gnpsdb, file = paste(input_dir,"gnps.rda", sep = ""))
        
        # delete the database in its format to free up space
        system(paste("rm", (paste(input_dir, "ALL_GNPS.mgf", sep = '')), sep = " "))
        
        writeLines(paste("GNPS saved at", Sys.time(), sep=" "),con=file.conn)
        
    }
    
    #mbank
    if (db == "all" || db =="mbank"){
        
        print("MassBank WORKS")
        
        page <- read_html("https://github.com/MassBank/MassBank-data/releases")
        page %>%
            html_nodes("a") %>%       # find all links
            html_attr("href") %>%     # get the url
            str_subset("MassBank_NIST.msp") -> tmp # find those that have the name MassBank_NIST.msp
        
        #download file
        system(paste("wget ",
                     "https://github.com", tmp[1], 
                     sep =  ""))
        
        mbank <- Spectra(paste(input_dir, "MassBank_NIST.msp", sep = ''), source = MsBackendMsp())
        save(mbank, file = paste(input_dir,"mbankNIST.rda", sep = ""))
        
        # delete the database in its format to free up space
        system(paste("rm", (paste(input_dir, "MassBank_NIST.msp", sep = '')), sep = " "))
        
        # obtain the month and year for the database release to add to summary
        res <- str_match(tmp[1], "download/\\s*(.*?)\\s*/MassBank_NIST")
        
        writeLines(paste("MassBank saved at", Sys.time(), "with release version", res[,2], sep=" "),con=file.conn)
    }
    
    # hmdb
    if (db == "all" || db =="hmdb"){
        
        print("HMDB WORKS")
        
        
        
        ####### Version Control ######
        
        # extract HMDB Current version
        html <- read_html("https://hmdb.ca/downloads")
        strings <- html%>% html_elements("a") %>% html_text2()
        ls <- unique(strings)
        hmdb_curr_ver <- c()
        for (i in ls){
            if (grepl("Current", i)){
            hmdb_curr_ver<- c(i, hmdb_curr_ver)
            }
        }
        
        
        
        
        ####### Download and unzip ######
        
        #Download file predicted MSMS spectra
        system(paste("wget",
                     "https://hmdb.ca/system/downloads/current/spectral_data/spectra_xml/hmdb_predicted_msms_spectra.zip",
                     sep = " "))
        # unzip
        system(paste("unzip", "hmdb_predicted_msms_spectra.zip", "-d",  paste(input_dir, "hmdb_predicted_msms_spectra", sep = ""), sep = " "))
    
        
        #Download file experimental MSMS spectra
        system(paste("wget",
                     "https://hmdb.ca/system/downloads/current/spectral_data/spectra_xml/hmdb_experimental_msms_spectra.zip",
                     sep = " "))
        # unzip
        system(paste("unzip", "hmdb_experimental_msms_spectra.zip", "-d", paste(input_dir, "hmdb_experimental_msms_spectra", sep = ""), sep = " "))
        
        
        
        
        ####### Load spectra in MsBackend #######
        
        hmdb_predfiles <- list.files(path = paste(input_dir, "hmdb_predicted_msms_spectra", sep = ''), full.names = TRUE)
        
        hmdb_predicted <- c()
        for (i in hmdb_predfiles){
            # load the spectra into MsBackendHMDB
            hmdb_pred <- Spectra(i, source = MsBackendHmdbXml())
            hmdb_predicted <- c(hmdb_predicted, hmdb_pred)
        }
        
        hmdb_expfiles <- list.files(path = paste(input_dir, "hmdb_experimental_msms_spectra", sep = ''), full.names = TRUE)
        
        hmdb_experimental <- c()
        
        for (j in hmdb_expfiles){
            # load the spectra into MsBackendHMDB
            hmdb_exp <- Spectra(j, source = MsBackendHmdb())
            hmdb_experimental <- c(hmdb_experimental, hmdb_exp)
        }
        
        
        hmdb <- hmdb_predicted + hmdb_experimental
        save(hmdb, file = paste(input_dir,"hmdb.rda", sep = ""))
        
        
        
        
        
        ####### Remove the XML files #######
        
        # delete the database in its format to free up space
        system(paste("rm -r", (paste(input_dir, "hmdb_predicted_msms_spectra", sep = '')), sep = " "))
        system(paste("rm -r", (paste(input_dir, "hmdb_experimental_msms_spectra", sep = '')), sep = " "))
        
        
        writeLines(paste("HMDB saved at", Sys.time(), "with release version", hmdb_curr_ver, sep=" "),con=file.conn)
    }
    
    #wrong input error message
    else if (!grepl(db, databases, fixed = TRUE)){
        stop("Wrong db input. Following inputs apply: gnps, hmdb, mbank or all")
    }
    close(file.conn)
    end_time <- Sys.time()
    print(end_time - start_time)
}


download_specDB(input_dir, db, error = TRUE)



