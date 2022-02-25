#! /usr/bin/Rscript

#' @title Module: Spectral DB Dereplication

#'
#' @description 

#' This module performs preprocessing on MS2 spectra from input 
#' .mzML files and from the spectra objects stored in spectral DBs 
#' downloaded as .rda object using the function download_specDB. The 
#' module then performs spectral database dereplication.
#' This module consists of spec2_Processing, peakdf, low_int, norm_int,
#' Listed below are the parameters for spec_dereplication

#' @param x precursor m/z
#'
#' @param db is either one of the spectral libraries which can be 
#'        gnps, hmdb, mbank or all
#'
#' @param result_dir generated from function ms2_rfilename for that result directory
#'
#' @param file_id generated from function ms2_rfilename for that result directory
#'
#' @param input_dir is input directory with all input files
#'
#' @param ppmx is the allowed ppm range between the m/z of the two spectra
#' being compared, by default it is 15

#' @return
#' Each best match generates a mirrorspectra of the match in the following created
#' 
#' A final csv that is a dataframe of results
#' 
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' # spec_dereplication(425.90, "all", "usr/project/file1.mzML", "id_1", "/usr/project/", ppmx =15)



# ---------- Preparations ----------
# Load libraries
library("Spectra")
library("MsBackendMgf")
library("MsBackendMsp")
library("MsBackendHmdb")
library("stringr")
library("MsCoreUtils")
library("readr")
library("dplyr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

x <- as.numeric(args[1])
db <- as.character(args[2])
result_dir <- as.character(args[3])
file_id <- as.character(args[4])
input_dir <- as.character(args[5])
ppmx <- as.numeric(args[6])
error <- as.logical(args[7])


# ---------- spec_dereplication ----------

# x is one pre_mz, db is GNPS, HMDB, MassBank
spec_dereplication<- function(x, db, result_dir, file_id, input_dir, ppmx, error = TRUE){
    
    ####-------------------------------------------------------------
    #### Define dependent Functions ----
    ####-------------------------------------------------------------
    
    #' Define a function to *normalize* the intensities
    norm_int <- function(y, ...) {
        maxint <- max(y[, "intensity"], na.rm = TRUE)
        y[, "intensity"] <- 100 * y[, "intensity"] / maxint
        y
    }
    #' Define a filtering function and remove peaks less than 0.05 of intensity
    low_int <- function(c, ...) {
        c > max(c, na.rm = TRUE) * 0.05
    }
    #' Specifying a function to draw peak labels
    label_fun <- function(p) {
        ints <- unlist(intensity(p))
        mzs <- format(unlist(mz(p)), digits = 4)
        mzs[ints < 5] <- ""
        mzs
    }
    
    
    
    #' processing on spectra with one precursor mass
    # inputs: 
    # x is precursor mass, 
    # spec is the spectra file (sps_all is mzML input processed spectra, gnps, hmdb or mbank), 
    # ppmx is ppm value
    spec2_Processing <- function(z, spec = "sps_all", ppmx = 15){
        if (spec == "sps_all"){
            #' Subset the dataset to MS2 spectra matching the m/z
            sps <- filterPrecursorMz(sps_all, mz = z + ppm(c(-z, z), 10))
        } else if (spec == "gnps"){
            #gnps spectra that contains precursor mass
            has_mz <- containsMz(gnps, mz = z, ppm = ppmx)
            #' Subset the GNPS Spectra
            sps <- gnps[has_mz]
        } else if (spec == "hmdb"){
            #hmdb spectra that contains precursor mass
            has_mz <- containsMz(hmdb, mz = z, ppm = ppmx)
            #' Subset the HMDB Spectra
            sps <- hmdb[has_mz]
        } else if (spec == "mbank"){
            has_mz <- containsMz(mbank,mz = z, ppm = ppmx)
            #' Subset the HMDB Spectra
            sps <- mbank[has_mz]
        }
    
        #wrong input error message
        else if (!grepl(db, databases, fixed = TRUE)){
            stop("Wrong db input. Following inputs apply: gnps, hmdb, mbank or all")
        }

        if (length(sps)>0){
            #' Apply the function to filter the spectra
            sps <- filterIntensity(sps, intensity = low_int)
            #' *Apply* the function to the data
            sps <- addProcessing(sps, norm_int)
            # cleaning peaks that are heavier or equal to the precursor mass
            pkd <- peaksData(sps)@listData
        
            #obtain the list of peaks that are higher or equal to precursor mass
            # y is peaksData from spectra
            removePrecursorPeaks <- function(y){
                y <- y[y[, "mz"] < z, ]
            }
            # use lapply to apply the function to the list of peaksData
            pkd <- lapply(pkd, removePrecursorPeaks)
            #store the indices of spectra with 0 peaks
            store_i <- c()
            for (i in 1:length(pkd)){
                if (is.null(nrow(pkd[[i]]))){
                    #convert the object of one peak into a matrix
                    mz <- pkd[[i]][[1]]
                    intensity <- pkd[[i]][[2]]
                    mat <- cbind(mz, intensity)
                    pkd[[i]] <- mat
                }else if(nrow(pkd[[i]])==0){
                    #store indices with 0 peaks
                    store_i <- c(store_i, i)
                }
            }
            # if 0 peaks, remove the relevant spectra from gnps or hmdb or mbank
            if (!(is.null(store_i))){
                pkd <- pkd[-(store_i)]
                sps <- sps[-(store_i)]
                peaksData(sps@backend)<- pkd
            }else {
                peaksData(sps@backend)<- pkd
            }
            return(sps)
        }
        else {
            sps <- NULL
            return(sps)
        }
    }

    
    
    #' obtain peaksData for each spectral matching between query and database spectra
    #inputs a is best match from Database, b is best match from query spectra
    peakdf <- function(a, b, ppmx){

        #' obtain peaklists for both query spectra and best matched spectra from MassBank
        z<- peaksData(a)[[1]] #from GNPS/HMDB
        y <- peaksData(b)[[1]] #from query
        if (!(nrow(z)==0)){
        #' Since we used 15 ppm, so to find the range, calculate the mass error range
        range <-  (ppmx/1000000)*y[,"mz"]
        y <- cbind(y, range)
        low_range <- y[,"mz"]-y[,"range"] # low range of m/z
        high_range <- y[,"mz"]+y[,"range"] # high range of m/z
            y <- cbind(y, low_range, high_range)
        #from GNPS/HMDB/MassBank spectra
        mz.z <- c()
        intensity.z <- c()
        #from query spectra
        mz.y <- c()
        intensity.y <- c()
        #difference between their intensity
        diff <- c()
            
        #' for all rows of y
        for (m in 1:nrow(y)){
            #' for all rows of z
                
            for(j in 1:nrow(z)){
                    
                ###################################################################
        
                ## IFELSE Statement no.2 -- LOOP 1.1.1.1
                    
                #' if the m/z of MB Spectra is within the 20 ppm range, save difference between intensities
                if (y[m,"low_range"] <= z[j, "mz"] && z[j, "mz"] <= y[m,"high_range"]){
            
                    #GNPS/HMDB
                    mz_z <- as.numeric(z[j, "mz"])
                    mz.z <- c(mz.z, mz_z)
                        
                    intensity_z <- as.numeric(z[j, "intensity"])
                    intensity.z <- c(intensity.z, intensity_z)
                        
                    #QUERY
                    mz_y <- as.numeric(y[m, "mz"])
                    mz.y <- c(mz.y, mz_y)
                   
                    intensity_y <- as.numeric(y[m, "intensity"])
                    intensity.y <- c(intensity.y, intensity_y)
                    
                    #Difference between intensities
                    difference <- as.numeric(abs(z[j, "intensity"]-y[m, "intensity"]))
                    diff <- c(diff, difference)
                }
            }
        }
        df_peaklists <- cbind(mz.y, intensity.y, mz.z, intensity.z, diff)
        return(df_peaklists)
        }
        else{
            df_peaklists <- NULL
            return(df_peaklists)
        }
        #output is a dataframe with mz and intensity from db spectra and query spectra and their difference
    }
    ####-------------------------------------------------------------
    #### Dereplication with all or GNPS ----
    ####-------------------------------------------------------------
    
    databases <- 'gnps, hmdb, mbank, all'
    
    if (db == "all" || db =="gnps"){
        nx <- 0
        nx <- nx+1
       
        #### input spec with pre_mz
        sps <- spec2_Processing(x, spec = "sps_all", ppmx)
        
        #### GNPS spec with pre_mz
        gnps_with_mz <- spec2_Processing(x, spec = "gnps", ppmx)
        
        dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/GNPS/", sep = ""), "./"), sep ="")
        if (!file.exists(dir_name)){
            dir.create(dir_name, recursive = TRUE)
        }
        
        if (length(sps) > 1 && length(gnps_with_mz) >1){
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
            
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            gnps_best_match <- gnps_with_mz[idx[2]]
            df_peaklists <- peakdf(gnps_best_match, sps[idx[1]], ppmx)
            
            if (!(is.null(df_peaklists))){
                
                #print("more spectra and more gnps spectra")
                
                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", gnps_best_match$SPECTRUMID, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps[idx[1]], gnps_with_mz[idx[2]], tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- max(res)
                GNPSmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                GNPSintScore <- mean(1-(df_peaklists[,"diff"]/100))
                GQMatchingPeaks <- nrow(df_peaklists)
                GNPSTotalPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                gQueryTotalPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                GNPSSMILES <- gnps_best_match$SMILES
                GNPSspectrumID <- gnps_best_match$SPECTRUMID
                GNPScompound_name <- gnps_best_match$NAME
                GNPSmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "GNPS"
            }
            else{
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- NA
                GNPSmzScore <- NA
                GNPSintScore <- NA
                GQMatchingPeaks <- NA
                GNPSTotalPeaks <- NA
                gQueryTotalPeaks<- NA
                GNPSSMILES <- NA
                GNPSspectrumID <- NA
                GNPScompound_name <- NA
                GNPSmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(gnps_with_mz) >1){
            
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            
            gx <- which(res == max(res))
            gx <- gx[1]
            gnps_best_match <- gnps_with_mz[gx]

            df_peaklists <- peakdf(gnps_best_match, sps, ppmx)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", gnps_best_match$SPECTRUMID, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, gnps_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off() 
                
                #' Identify the best-matching pair
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- max(res)
                GNPSintScore <- mean(1-(df_peaklists[,"diff"]/100))
                GNPSmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                GQMatchingPeaks <- nrow(df_peaklists)
                GNPSTotalPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                gQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                GNPScompound_name <- gnps_best_match$NAME 
                GNPSspectrumID <- gnps_best_match$SPECTRUMID
                GNPSSMILES <- gnps_best_match$SMILES
                GNPSmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "GNPS"
            }
            else{
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- NA
                GNPSmzScore <- NA
                GNPSintScore <- NA
                GQMatchingPeaks <- NA
                GNPSTotalPeaks <- NA
                gQueryTotalPeaks<- NA
                GNPSSMILES <- NA
                GNPSspectrumID <- NA
                GNPScompound_name <- NA
                GNPSmirrorSpec <- NA
                Source <- NA
            }
        
        }
        else if (length(sps) > 1 && length(gnps_with_mz) == 1){
            
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
            #' obtain MB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(gnps_with_mz, sps[gx], ppmx)
            gnps_best_match <- gnps_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", gnps_best_match$SPECTRUMID, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, gnps_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- max(res)
                GNPSintScore <- mean(1-(df_peaklists[,"diff"]/100))
                GNPSmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                GQMatchingPeaks <- nrow(df_peaklists)
                GNPSTotalPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                gQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                GNPScompound_name <- gnps_best_match$NAME 
                GNPSspectrumID <- gnps_best_match$SPECTRUMID
                GNPSSMILES <- gnps_best_match$SMILES
                GNPSmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "GNPS"
                }
            else{
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- NA
                GNPSmzScore <- NA
                GNPSintScore <- NA
                GQMatchingPeaks <- NA
                GNPSTotalPeaks <- NA
                gQueryTotalPeaks<- NA
                GNPSSMILES <- NA
                GNPSspectrumID <- NA
                GNPScompound_name <- NA
                GNPSmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(gnps_with_mz) == 1){
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
            gnps_best_match <- gnps_with_mz
            df_peaklists <- peakdf(gnps_best_match, sps, ppmx)
            if (!(is.null(df_peaklists))){
                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", gnps_best_match$SPECTRUMID, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, gnps_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- max(res)
                GNPSintScore <- mean(1-(df_peaklists[,"diff"]/100))
                GNPSmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                GQMatchingPeaks <- nrow(df_peaklists)
                GNPSTotalPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                gQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                GNPScompound_name <- gnps_best_match$NAME 
                GNPSspectrumID <- gnps_best_match$SPECTRUMID
                GNPSSMILES <- gnps_best_match$SMILES
                GNPSmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "GNPS"
            }
            else{
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                GNPSmax_similarity <- NA
                GNPSmzScore <- NA
                GNPSintScore <- NA
                GQMatchingPeaks <- NA
                GNPSTotalPeaks <- NA
                gQueryTotalPeaks<- NA
                GNPSSMILES <- NA
                GNPSspectrumID <- NA
                GNPScompound_name <- NA
                GNPSmirrorSpec <- NA
                Source <- NA
            }
                
        }
        else{
            id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
            premz <- x
            rtmin <- max(sps$rtime)
            rtmax <- min(sps$rtime)
            rtmed <- median(sps$rtime, na.rm = TRUE)
            rtmean <- mean(sps$rtime, na.rm = TRUE)
            GNPSmax_similarity <- NA
            GNPSmzScore <- NA
            GNPSintScore <- NA
            GQMatchingPeaks <- NA
            GNPSTotalPeaks <- NA
            gQueryTotalPeaks<- NA
            GNPSSMILES <- NA
            GNPSspectrumID <- NA
            GNPScompound_name <- NA
            GNPSmirrorSpec <- NA
            Source <- NA

        }
        
        df_gnps <- cbind(id_X, premz, rtmin, rtmax, rtmed, rtmean, GNPSmax_similarity, GNPSmzScore, 
                        GNPSintScore, GQMatchingPeaks, GNPSTotalPeaks, gQueryTotalPeaks, 
                        GNPSSMILES, GNPSspectrumID, GNPScompound_name, GNPSmirrorSpec, Source)
        
        #write.csv(df_gnps, str_remove(paste(input_dir, result_dir, "/spectral_dereplication/gnps.csv", sep = ""), "./"))
        write.csv(df_gnps, paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/gnps.csv", sep = ""), "./"), sep = ""))
        #return(df_gnps)
        
    }
    
    ####-------------------------------------------------------------
    #### Dereplication with all or HMDB ----
    ####-------------------------------------------------------------
    
    if (db == "all" || db =="hmdb"){
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(x, spec = "sps_all", ppmx = NULL)
        
        #### HMDB spec with pre_mz
        hmdb_with_mz <- spec2_Processing(x, spec = "hmdb", ppmx = 15)
        
        
        dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/HMDB/", sep = ""), "./"), sep ="")
        if (!file.exists(dir_name)){
            dir.create(dir_name, recursive = TRUE)
        }
        if (length(sps) > 1 && length(hmdb_with_mz) > 1){
            #' Compare experimental spectra against HMDB
            res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)
            
            #' obtain HMDB spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            hmdb_best_match <- hmdb_with_mz[idx[2]]
            df_peaklists <- peakdf(hmdb_best_match, sps[idx[1]], ppmx)
            
            if (!(is.null(df_peaklists))){

                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", hmdb_best_match$compound_id, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps[idx[1]], hmdb_with_mz[idx[2]], tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                HMDBcompoundID <- hmdb_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
            }
            else{
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- "none"
                HMDBmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(hmdb_with_mz) >1){
            #' Compare experimental spectra against HMDB
            res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)
            
            #' obtain HMDB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            hmdb_best_match <- hmdb_with_mz[gx]

            df_peaklists <- peakdf(hmdb_best_match, sps, ppmx)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", hmdb_best_match$compound_id, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, hmdb_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off() 
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                HMDBcompoundID <- hmdb_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
            }
            else{
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- "none"
                HMDBmirrorSpec <- NA
                Source <- NA
            }
        
        }
        else if (length(sps) > 1 && length(hmdb_with_mz) == 1){
            #' Compare experimental spectra against HMDB
            res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)
            #' obtain hmdb spectra that matches the most with m/z MS2 spectra
            
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(hmdb_with_mz, sps[gx], ppmx)
            hmdb_best_match <- hmdb_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){

                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", hmdb_best_match$compound_id, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, hmdb_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                HMDBcompoundID <- hmdb_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
                }
            else{

                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- "none"
                HMDBmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(hmdb_with_mz) == 1){
            #' Compare experimental spectra against HMDB
            res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)
            hmdb_best_match <- hmdb_with_mz
            df_peaklists <- peakdf(hmdb_best_match, sps, ppmx)
            if (!(is.null(df_peaklists))){
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", hmdb_best_match$compound_id, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, hmdb_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                HMDBcompoundID <- hmdb_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
            }
            else{

                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- "none"
                HMDBmirrorSpec <- NA
                Source <- NA
            }
                
        }
        else{

            id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
            premz <- x
            rtmin <- max(sps$rtime)
            rtmax <- min(sps$rtime)
            rtmed <- median(sps$rtime, na.rm = TRUE)
            rtmean <- mean(sps$rtime, na.rm = TRUE)
            HMDBmax_similarity <- NA
            HMDBmzScore <- NA
            HMDBintScore <- NA
            HQMatchingPeaks <- NA
            HMDBTotalPeaks <- NA
            hQueryTotalPeaks<- NA
            HMDBcompoundID <- "none"
            HMDBmirrorSpec <- NA
            Source <- NA
        }
        
        df_hmdb <- cbind(id_X, premz, rtmin, rtmax, rtmed, rtmean, HMDBmax_similarity, HMDBmzScore, 
                         HMDBintScore, HQMatchingPeaks, HMDBTotalPeaks, hQueryTotalPeaks, 
                         HMDBcompoundID, HMDBmirrorSpec, Source)
        
        #write.csv(df_hmdb, str_remove(paste(result_dir, "/spectral_dereplication/hmdb.csv", sep = ""), "./"))
        write.csv(df_hmdb, paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/hmdb.csv", sep = ""), "./"), sep = ""))
        #return(df_hmdb)
    }
    
    ####-------------------------------------------------------------
    #### Dereplication with all or MassBank ----
    ####-------------------------------------------------------------
    
    if (db == "all" || db =="mbank"){
        
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(x, spec = "sps_all", ppmx = NULL)
        
        #### GNPS spec with pre_mz
        mbank_with_mz <- spec2_Processing(x, spec = "mbank", ppmx = 15)
        
        dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/MassBank/", sep = ""), "./"), sep ="")
        if (!file.exists(dir_name)){
            dir.create(dir_name, recursive = TRUE)
        }
        if (length(sps) > 1 && length(mbank_with_mz) >1){
            
            #' Compare experimental spectra against MassBank
            res <- compareSpectra(sps, mbank_with_mz, ppm = 15)
            
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            mbank_best_match <- mbank_with_mz[idx[2]]
            df_peaklists <- peakdf(mbank_best_match, sps[idx[1]], ppmx)
            
            if (!(is.null(df_peaklists))){
                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", mbank_best_match$accession, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps[idx[1]], mbank_with_mz[idx[2]], tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- max(res)
                MBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                MBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                MQMatchingPeaks <- nrow(df_peaklists)
                MBTotalPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                mQueryTotalPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                MBformula <- mbank_best_match$formula
                MBinchiKEY <- mbank_best_match$inchikey
                MBspectrumID <- mbank_best_match$accession
                MBcompound_name <- mbank_best_match$name
                MBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "MassBank"
            }
            else{

                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- NA
                MBmzScore <- NA
                MBintScore <- NA
                MQMatchingPeaks <- NA
                MBTotalPeaks <- NA
                mQueryTotalPeaks<- NA
                MBformula <- NA
                MBinchiKEY <- NA
                MBspectrumID <- NA
                MBcompound_name <- NA
                MBmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(mbank_with_mz) >1){
            
            #' Compare experimental spectra against MassBank
            res <- compareSpectra(sps, mbank_with_mz, ppm = 15)
            #' obtain MassBank spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            mbank_best_match <- mbank_with_mz[gx]

            df_peaklists <- peakdf(mbank_best_match, sps, ppmx)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){

                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", mbank_best_match$accession, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, mbank_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off() 
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- max(res)
                MBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                MBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                MQMatchingPeaks <- nrow(df_peaklists)
                MBTotalPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                mQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                MBformula <- mbank_best_match$formula
                MBinchiKEY <- mbank_best_match$inchikey
                MBspectrumID <- mbank_best_match$accession
                MBcompound_name <- mbank_best_match$name
                MBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "MassBank"
            }
            else{

                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- NA
                MBmzScore <- NA
                MBintScore <- NA
                MQMatchingPeaks <- NA
                MBTotalPeaks <- NA
                mQueryTotalPeaks<- NA
                MBformula <- NA
                MBinchiKEY <- NA
                MBspectrumID <- NA
                MBcompound_name <- NA
                MBmirrorSpec <- NA
                Source <- NA
            }
        
        }
        else if (length(sps) > 1 && length(mbank_with_mz) == 1){
            
            #' Compare experimental spectra against MassBank
            res <- compareSpectra(sps, mbank_with_mz, ppm = 15)
            #' obtain MB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(mbank_with_mz, sps[gx], ppmx)
            mbank_best_match <- mbank_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){

                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", mbank_best_match$accession, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, mbank_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- max(res)
                MBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                MBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                MQMatchingPeaks <- nrow(df_peaklists)
                MBTotalPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                mQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                MBformula <- mbank_best_match$formula
                MBinchiKEY <- mbank_best_match$inchikey
                MBspectrumID <- mbank_best_match$accession
                MBcompound_name <- mbank_best_match$name
                MBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "MassBank"
                }
            else{

                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- NA
                MBmzScore <- NA
                MBintScore <- NA
                MQMatchingPeaks <- NA
                MBTotalPeaks <- NA
                mQueryTotalPeaks<- NA
                MBformula <- NA
                MBinchiKEY <- NA
                MBspectrumID <- NA
                MBcompound_name <- NA
                MBmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(mbank_with_mz) == 1){
            
            #' Compare experimental spectra against MassBank
            res <- compareSpectra(sps, mbank_with_mz, ppm = 15)
            mbank_best_match <- mbank_with_mz
            df_peaklists <- peakdf(mbank_best_match, sps, ppmx)
            if (!(is.null(df_peaklists))){

                
                #' plotMirror
                name_plotmirror <- paste(dir_name, x,"_spectra_vs_", mbank_best_match$accession, "_spectra.pdf", sep ="")
                pdf(name_plotmirror)
                plotSpectraMirror(sps, mbank_best_match, tolerance = 0.2,
                                    labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                    labelSrt = -30)
                grid()
                dev.off()
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- max(res)
                MBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                MBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                MQMatchingPeaks <- nrow(df_peaklists)
                MBTotalPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                mQueryTotalPeaks<- nrow(peaksData(sps)[[1]])
                MBformula <- mbank_best_match$formula
                MBinchiKEY <- mbank_best_match$inchikey
                MBspectrumID <- mbank_best_match$accession
                MBcompound_name <- mbank_best_match$name
                MBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "MassBank"
            }
            else{

                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                rtmean <- mean(sps$rtime, na.rm = TRUE)
                MBmax_similarity <- NA
                MBmzScore <- NA
                MBintScore <- NA
                MQMatchingPeaks <- NA
                MBTotalPeaks <- NA
                mQueryTotalPeaks<- NA
                MBformula <- NA
                MBinchiKEY <- NA
                MBspectrumID <- NA
                MBcompound_name <- NA
                MBmirrorSpec <- NA
                Source <- NA
            }
                
        }
        else{

            id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
            premz <- x
            rtmin <- max(sps$rtime)
            rtmax <- min(sps$rtime)
            rtmed <- median(sps$rtime, na.rm = TRUE)
            rtmean <- mean(sps$rtime, na.rm = TRUE)
            MBmax_similarity <- NA
            MBmzScore <- NA
            MBintScore <- NA
            MQMatchingPeaks <- NA
            MBTotalPeaks <- NA
            mQueryTotalPeaks<- NA
            MBformula <- NA
            MBinchiKEY <- NA
            MBspectrumID <- NA
            MBcompound_name <- NA
            MBmirrorSpec <- NA
            Source <- NA

        }
        
        df_mbank <- cbind(id_X, premz, rtmin, rtmax, rtmed, rtmean, MBmax_similarity, MBmzScore, MBintScore, MQMatchingPeaks, 
                         MBTotalPeaks, mQueryTotalPeaks, MBformula, MBinchiKEY, MBspectrumID, MBcompound_name, MBmirrorSpec, 
                         Source)
        write.csv(df_mbank, paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/mbank.csv", sep = ""), "./"), sep = ""))
        
        #return(df_mbank)
    }
    #wrong input error message
    else if (!grepl(db, databases, fixed = TRUE)){
        stop("Wrong db input. Following inputs apply: gnps, hmdb, mbank or all")
    }
    
}
spec_dereplication(x, db, result_dir, file_id, input_dir, ppmx)

