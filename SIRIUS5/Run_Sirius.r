# provenance library rdtLite integration
# library(rdtLite)
# options(prov.dir = "./prov_metfrag", snapshot.size = 10000)
# prov.init(prov.dir = ".")

#define arguments
args = commandArgs(trailingOnly = TRUE)

# Start time
start.time <- Sys.time()

files <- args[1]
profile <- args[2]
db <- args[3]

run_sirius <- function(files, ppm_max = 5, ppm_max_ms2 = 15, QC, SL , SL_path, candidates = 30, profile, db){
    files <- read.csv(files, sep = "\t")

    for (b in 1:nrow(files)){

        if (QC){
            if (is.na(files[b, "isotopes"])){
                system(paste("sirius --input", files[b, "sirius_param_file"],
                             "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,
                             "canopus write-summaries --output", files[b, "outputNames"],
                             sep = " "))
                if(!file.exists(files[b, "outputNames"])){
                    system(paste("sirius --input", files[b, "sirius_param_file"],
                                 "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,
                                 "canopus write-summaries --output", files[b, "outputNames"],
                                 sep = " "))
                }
                if(SL){
                    system(paste("sirius --input", files[b, "sirius_param_file"],
                             "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database",  SL_path ,
                             "canopus write-summaries --output", files[b, "outputNamesSL"],
                             sep = " "))
                }

            }
            else if(files[b, "isotopes"] == "present"){
                system(paste("sirius --input", files[b, "sirius_param_file"],
                             "formula --profile", profile, "--candidates", candidates, "--ppm-max", ppm_max,"--ppm-max-ms2", ppm_max_ms2,"fingerprint structure --database", db ,
                             "canopus write-summaries --output", files[b, "outputNames"],
                             sep = " "))
                if(!file.exists(files[b, "outputNames"])){
                     system(paste("sirius --input", files[b, "sirius_param_file"],
                                 "formula --profile", profile, "--candidates", candidates, "--ppm-max", ppm_max,"--ppm-max-ms2", ppm_max_ms2,"fingerprint structure --database", db ,
                                 "canopus write-summaries --output", files[b, "outputNames"],
                                 sep = " "))
                }
                if(SL){
                    system(paste("sirius --input", files[b, "sirius_param_file"],
                             "formula --profile", profile, "--candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database",  SL_path ,
                             "canopus write-summaries --output", files[b, "outputNamesSL"],
                             sep = " "))
                }
            }
        }
        else{
            system(paste("sirius --input", files[b, "sirius_param_file"], 
                        "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,
                        "canopus write-summaries --output", files[b, "outputNames"],
                        sep = " "))
            if(!file.exists(files[b, "outputNames"])){
                system(paste("sirius --input", files[b, "sirius_param_file"],
                             "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database", db ,
                             "canopus write-summaries --output", files[b, "outputNames"],
                             sep = " "))
            }
            if(SL){
                system(paste("sirius --input", files[b, "sirius_param_file"],
                            "formula --profile", profile, "--no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"fingerprint structure --database",  SL_path ,
                            "canopus write-summaries --output", files[b, "outputNamesSL"],
                            sep = " "))
            }
        }
        Sys.sleep(5)

    }


}

run_sirius(files, ppm_max = 5, ppm_max_ms2 = 15, QC = FALSE, SL = FALSE, SL_path = NA, candidates = 30, profile, db)

end.time <- Sys.time()

# Time taken to run the analysis for MAW-R
time.taken <- end.time - start.time
print(time.taken)