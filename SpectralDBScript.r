# Load Libraries
library(Spectra)
library(MsBackendMgf)
library(MsBackendHmdb)
library(MsCoreUtils)
library(MsBackendMsp)

# Load functions file
input_dir <- paste(getwd(), "/", sep = '')
load(file = paste(input_dir, "R_Functions.RData", sep = ''))
# Download the databases
download_specDB(input_dir, db = "all")


