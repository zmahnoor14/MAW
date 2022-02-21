#! /usr/bin/Rscript

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

input_dir <- as.character(args[1])
db <- as.character(args[2])
error <- as.logical(args[3])

download_specDB <- function(input_dir, db = "all", error = TRUE){
    # gnps
    if (db == "all" || db =="gnps"){
        print("gnps is ok")
    }
    # gnps
    if (db == "all" || db =="hmdb"){
        print("hmdb is ok")
    }
    # gnps
    if (db == "all" || db =="mbank"){
        print("mbank is ok")
    }
}
download_specDB(input_dir, db, error = TRUE)


