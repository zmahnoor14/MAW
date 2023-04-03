library(stringr)
library(readr)
library(dplyr)

start.time <- Sys.time()

# ---------- Script ----------
# input directory
input_dir <- paste(getwd(), "/data", sep = "")
input_dir

run_sirius <- function(files, ppm_max = 5, ppm_max_ms2 = 15, QC = TRUE, SL = TRUE, SL_path, candidates = 30, profile, db){
    files <- read.csv(files, sep = "\t")

    for (b in 1:nrow(files)){

        if (QC){
            if (is.na(files[b, "isotopes"])){
                system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                             "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,"canopus",
                             sep = " "))
                if(!file.exists(files[b, "outputNames"])){
                    system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                                 "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,"canopus",
                                 sep = " "))
                }
                if(SL){
                    system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNamesSL"],
                             "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database",  SL_path ,"canopus",
                             sep = " "))
                }

            }
            else if(files[b, "isotopes"] == "present"){
                system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                             "formula --profile", profile, "--candidates", candidates, "--ppm-max", ppm_max,"--ppm-max-ms2", ppm_max_ms2,"fingerprint structure --database", db ,"canopus",
                             sep = " "))
                if(!file.exists(files[b, "outputNames"])){
                     system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                                 "formula --profile", profile, "--candidates", candidates, "--ppm-max", ppm_max,"--ppm-max-ms2", ppm_max_ms2,"fingerprint structure --database", db ,"canopus",
                                 sep = " "))
                }
                if(SL){
                    system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNamesSL"],
                             "formula --profile", profile, "--candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database",  SL_path ,"canopus",
                             sep = " "))
                }
            }
        }
        else{
            system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                        "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,"canopus",
                        sep = " "))
            if(!file.exists(files[b, "outputNames"])){
                system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                             "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,"canopus",
                             sep = " "))
            }
            if(SL){
                system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNamesSL"],
                            "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database",  SL_path ,"canopus",
                            sep = " "))
            }
        }
        Sys.sleep(5)

    }

}






# please change parameters according to your experiment

run_sirius(files= './insilico/MS1DATA_SiriusP.tsv',
           ppm_max = 5, 
           ppm_max_ms2 = 15, 
           QC = FALSE, 
           SL = FALSE, 
           SL_path = NA, 
           candidates = 30, 
           profile = "orbitrap", 
           db = "coconut")


