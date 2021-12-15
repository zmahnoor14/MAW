library(Spectra)
library(MsBackendMgf)
library(MsBackendHmdb)
library(MsCoreUtils)
library(MsBackendMsp)
library(stringr)

input_dir <- paste(getwd(), "/", sep = '')

input_dir

# Load all spectral libraries
load(file = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/gnps.rda")
load(file = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/hmdb.rda")
load(file = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/mbankNIST.rda")

load(file = paste(input_dir, "R_Functions.RData", sep = ''))

input_table <- ms2_rfilename(input_dir)

input_table

x = "/Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_propanoyl_carnitine.mzML"

precursorMZs <- spec_Processing(x)
sps_all <- precursorMZs[[1]]
pre_mz<- precursorMZs[[2]]

a <- precursorMZs[[2]]

df_gnps <- spec_dereplication(a, 
                              db = "GNPS", 
                              result_dir = "/Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_propanoyl_carnitine",
                             file_id = "file_13",
                             input_dir)

df_hmdb <- spec_dereplication(a, 
                              db = "HMDB", 
                              result_dir = "/Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_propanoyl_carnitine",
                             file_id = "file_13",
                             input_dir)

df_hmdb

df_gnps

merge(df_hmdb, df_gnps)







sps <- spec2_Processing(a, spec = "sps_all", ppmx = NULL)

sps

gnps_with_mz <- spec2_Processing(a, spec = "gnps", ppmx = 15)

res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)

max(res)

idx <- which(res == max(res), arr.ind = TRUE)

gnps_best_match <- gnps_with_mz[idx[2]]

spectraVariables(gnps_best_match)

gnps_best_match$NAME


