#! /usr/local/bin/Rscript
args <- commandArgs(trailingOnly=TRUE)
print(args)
input_dir <- as.character(args[1])
db <- as.character(args[2])

say_Hello <- function(input_dir, db){
    if (db == "gnps"){
        system(paste("wget -P", 
                 input_dir,
                 "https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.mgf",
                 sep =  " "))
    }
}
say_Hello(input_dir, db)


