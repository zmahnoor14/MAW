library(Spectra)
library(MsBackendMgf)
library(MsBackendHmdb)
library(MsCoreUtils)
library(MsBackendMsp)
library(stringr)
library(readr)
library(dplyr)

input_dir <- paste(getwd(), "/", sep = '')
input_dir

# Load all spectral libraries
load(file = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/gnps.rda")
load(file = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/hmdb.rda")
load(file = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/mbankNIST.rda")

load(file = paste(input_dir, "R_Functions.RData", sep = ''))

input_table <- ms2_rfilename(input_dir)

input_table

load(file = paste(input_dir, "R_Functions.RData", sep = ''))

precursorMZs <- spec_Processing(as.character(input_table[3, "mzml_files"]))
sps_all <- precursorMZs[[1]]
pre_mz<- precursorMZs[[2]]
for (a in pre_mz){
    df_mbank <- spec_dereplication(a, 
                            db = "MassBank", 
                            result_dir = "/Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_methionine_sulfoxide",
                            file_id = "file_11",
                            input_dir)
    df_gnps <- spec_dereplication(a, 
                            db = "GNPS", 
                            result_dir = "/Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_methionine_sulfoxide",
                            file_id = "file_11",
                            input_dir)
    df_hmdb <- spec_dereplication(a, 
                            db = "HMDB", 
                            result_dir = "/Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_methionine_sulfoxide",
                            file_id = "file_11",
                            input_dir)
    }

df_mbank

load(file = paste(input_dir, "R_Functions.RData", sep = ''))

spec_pr <- spec_Processing(input_table[1, "mzml_files"])

spec_pr2 <- ms2_peaks(spec_pr, input_table[1, "ResultFileNames"])

spec_pr3 <- ms1_peaks(spec_pr2, y = NA, result_dir = input_table[1, "ResultFileNames"], QCfile = FALSE)

spec_pr3

sirius_param_files <- sirius_param(spec_pr3, result_dir = input_table[1, "ResultFileNames"])

sirius_param_files

system("sirius --input /Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_acetyl_carnitine/insilico/SIRIUS/1_NA_iso_NA_MS1p_204.122756958008_SIRIUS_param.ms --output /Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_acetyl_carnitine/insilico/SIRIUS/1_NA_iso_NA_MS1p_204.122756958008_SIRIUS_param.json formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates 30 --ppm-max 5 --ppm-max-ms2 15 structure --database ALL canopus")



sirius_postprocess(input_table[1, "ResultFileNames"])

sirius_postprocess(input_table[1, "ResultFileNames"], SL = FALSE)



AdductsMF <- read.csv(paste("/MetFrag_AdductTypes.csv", sep = ''))
