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

adducts = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetFrag_AdductTypes.csv"



input_table <- ms2_rfilename(input_dir)

input_table

for (i in nrow(input_table)){
    spec_pr <- spec_Processing(as.character(input_table[i, "mzml_files"]))
    sps_all <- spec_pr[[1]]
    pre_mz<- spec_pr[[2]]
    
    for (a in pre_mz){
        
        df_mbank <- spec_dereplication(a, 
                                db = "MassBank", 
                                result_dir = input_table[i, "ResultFileNames"],
                                file_id = input_table[i, "File_id"],
                                input_dir)
        df_gnps <- spec_dereplication(a, 
                                db = "GNPS", 
                                result_dir = input_table[i, "ResultFileNames"],
                                file_id = input_table[i, "File_id"],
                                input_dir)
        df_hmdb <- spec_dereplication(a, 
                                db = "HMDB", 
                                result_dir = input_table[i, "ResultFileNames"],
                                file_id = input_table[i, "File_id"],
                                input_dir)
    }
    
    spec_pr2 <- ms2_peaks(spec_pr, input_table[i, "ResultFileNames"])
    spec_pr3 <- ms1_peaks(spec_pr2, y = NA, result_dir = input_table[i, "ResultFileNames"], QCfile = FALSE)
    sirius_pfiles <- sirius_param(spec_pr3, result_dir = input_table[i, "ResultFileNames"], SL = FALSE)
    
    sirius_pfiles <- data.frame(sirius_pfiles)
    
    for (b in nrow(sirius_pfiles)){
        system(paste("sirius --input", sirius_pfiles[b, "sirius_param_file"], 
                     "--output", sirius_pfiles[b, "outputNames"],
                    "formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates 30 --ppm-max 5 --ppm-max-ms2 15 structure --database ALL canopus", 
                     sep = " "))
    }
    
    sirius_pproc <- sirius_postprocess(input_table[i, "ResultFileNames"], SL = FALSE)
    
    met_param <- metfrag_param(sirius_pproc, result_dir = input_table[i, "ResultFileNames"],
                               input_dir, adducts, sl_mtfrag = NA, SL = FALSE)
    for (files in met_param){
        system(paste("java -jar",  "/Users/mahnoorzulfiqar/MetFrag2.4.5-CL.jar", files))
    }
}

spec_pr <- spec_Processing(as.character(input_table[1, "mzml_files"]))
sps_all <- spec_pr[[1]]
pre_mz<- spec_pr[[2]]
    

spec_pr2 <- ms2_peaks(spec_pr, input_table[1, "ResultFileNames"])
spec_pr3 <- ms1_peaks(spec_pr2, y = NA, result_dir = input_table[1, "ResultFileNames"], QCfile = FALSE)
sirius_pfiles <- sirius_param(spec_pr3, result_dir = input_table[1, "ResultFileNames"], SL = FALSE)
    
sirius_pfiles <- data.frame(sirius_pfiles)
for (b in nrow(sirius_pfiles)){
    system(paste("sirius --input", sirius_pfiles[b, "sirius_param_file"], 
                    "--output", sirius_pfiles[b, "outputNames"],
                "formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates 30 --ppm-max 5 --ppm-max-ms2 15 structure --database ALL canopus", 
                    sep = " "))
}
    
sirius_pproc <- sirius_postprocess(input_table[1, "ResultFileNames"], SL = FALSE)
    
met_param <- metfrag_param(sirius_pproc, result_dir = input_table[1, "ResultFileNames"],
                            input_dir, adducts, sl_mtfrag = NA, SL = FALSE)
for (files in met_param){
    system(paste("java -jar",  "/Users/mahnoorzulfiqar/MetFrag2.4.5-CL.jar", files))
}

#system('java -jar /Users/mahnoorzulfiqar/MetFrag2.4.5-CL.jar /Users/mahnoorzulfiqar/Standards_CodeSet/VN_211016_acetyl_carnitine/insilico/MetFrag/1_id_M204R149ID1_mz_204.122756958008_rt_148.997391_db_PubChem.txt')






