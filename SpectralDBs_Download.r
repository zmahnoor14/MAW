# 4 dependencies for 
library(Spectra)
library(MsBackendMgf)
library(MsBackendHmdb)
library(MsBackendMsp)
# 3 dependencies for latest MassBank version
library(rvest)
library(stringr)
library(xml2)

## input directory ##
input_dir <- paste(getwd(), "/", sep = '')
input_dir
#use variable "input_dir"

download_specDB(input_dir, db = "all")

## Define function to download one of the open source spectral databases

## inputs
# input_dir = full directory where all MZML input files
# db = either one of the spectral libraries which can be gnps, hmdb, mbank or all

download_specDB <- function(input_dir, db){
    
    summaryFile <- paste(input_dir, "summaryFile.txt", sep = "")
    file.create(summaryFile, recursive = TRUE)
    file.conn <- file(summaryFile)
    open(file.conn, open = "at")
            
    
    if (db == "all" || db =="gnps"){
        # gnps
        
        # Download file
        system(paste("wget -P", 
                     input_dir,
                     "https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.mgf", 
                     sep =  " "))
        # load the spectra into MsBackendMgf
        gnps <- Spectra(paste(input_dir, "ALL_GNPS.mgf", sep = ''), source = MsBackendMgf())
        save(gnps, file = paste(input_dir,"gnps.rda", sep = ""))
        
        # delete the database in its format to free up space
        system(paste("rm", (paste(input_dir, "ALL_GNPS.mgf", sep = '')), sep = " "))
        
        writeLines(paste("GNPS saved at", Sys.time(), sep=" "),con=file.conn)
        
    }
    if (db == "all" || db =="hmdb"){
        # hmdb
        
        #Download file
        system(paste("wget - P", input_dir,
                     "https://hmdb.ca/system/downloads/current/spectral_data/spectra_xml/hmdb_all_spectra.zip",
                     sep = " "))
        # unzip
        system(paste("unzip -d", input_dir, paste(input_dir, "hmdb_all_spectra.zip", sep = ""), sep = " "))
        # load the spectra into MsBackendHMDB
        hmdb <- Spectra(paste(input_dir, "hmdb_all_spectra.xml", sep = ''), source = MsBackendHmdb())
        save(hmdb, file = paste(input_dir,"hmdb.rda", sep = ""))
        
        # delete the database in its format to free up space
        system(paste("rm", (paste(input_dir, "hmdb_all_spectra.xml", sep = '')), sep = " "))
        
        writeLines(paste("HMDB saved at", Sys.time(), sep=" "),con=file.conn)
    }
    if (db == "all" || db =="mbank"){
        #mbank
        
        page <- read_html("https://github.com/MassBank/MassBank-data/releases")
        page %>%
            html_nodes("a") %>%       # find all links
            html_attr("href") %>%     # get the url
            str_subset("MassBank_NIST.msp") -> tmp # find those that have the name MassBank_NIST.msp
        
        #download file
        system(paste("wget -P", input_dir,
                     "https://github.com/", tmp[1], 
                     sep =  " "))
        
        mbank <- Spectra(paste(input_dir, "MassBank_NIST.msp", sep = ''), source = MsBackendMgf())
        save(mbank, file = paste(input_dir,"mbankNIST.rda", sep = ""))
        
        # delete the database in its format to free up space
        system(paste("rm", (paste(input_dir, "MassBank_NIST.msp", sep = '')), sep = " "))
        
        # obtain the month and year for the database release to add to summary
        res <- str_match(tmp[1], "download/\\s*(.*?)\\s*/MassBank_NIST")
        
        writeLines(paste("MassBank saved at", Sys.time(), "with release version", res[,2], sep=" "),con=file.conn)
    }
    close(file.conn)
}

## outputs
# Spectral DB saved with the following name in the input_dir: gnps.rda, hmdb.rda, mbankNIST.rda
# summary file saved as summaryFile.txt which contains the timings of saved databases

## Usage:
# download_specDB(input_dir, db = "all")

load_specDB <- function(input_dir, db){
    if (db == "all" || db =="gnps"){
        load(file = paste(input_dir,"gnps.rda", sep = ""))
    }
    if (db == "all" || db =="hmdb"){
        load(file = paste(input_dir,"hmdb.rda", sep = ""))
    }
    if (db == "all" || db =="mbank"){
        load(file = paste(input_dir,"mbank.rda", sep = ""))
    }
}
