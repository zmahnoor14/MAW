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

#' @param pre_tbl txt file of list of precursor m/z (from func spec_Processing)
#'
#' @param proc_mzml processed mzML spectra file (from func spec_Processing)
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
#' A final csv that is a dataframe of results (gnps.csv, mbank.csv, hmdb.csv)
#' 
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' # spec_dereplication("/usr/project/file1/premz_list.txt", "/usr/project/file1/processedSpectra.mzmL","all", "/usr/project/file1.mzML", "id_1", "/usr/project/", ppmx =15)



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

pre_tbl <- as.character(args[1])
proc_mzml <- as.character(args[2])
db <- as.character(args[3])
result_dir <- as.character(args[4])
file_id <- as.character(args[5])
input_dir <- as.character(args[6])
ppmx <- as.numeric(args[7])
error <- as.logical(args[8])

# Functions spec_dereplication is dependent on


##-----------------------------------------------------------------
## filter intensity 
##-----------------------------------------------------------------

#' Define a filtering function and remove peaks less than 0.05 of intensity
low_int <- function(c, ...) {
    c > max(c, na.rm = TRUE) * 0.05
}

# Usage:
# filterIntensity(spectra_object, intensity = low_int)

#' filterIntensity is a predefined function in Spectra package

##-----------------------------------------------------------------
## normalize intensity 
##-----------------------------------------------------------------

#' Define a function to *normalize* the intensities
norm_int <- function(y, ...) {
    maxint <- max(y[, "intensity"], na.rm = TRUE)
    y[, "intensity"] <- 100 * y[, "intensity"] / maxint
    y
}




##-----------------------------------------------------------------
## Plotting Mirror Spectra 
##-----------------------------------------------------------------

#' Specifying a function to draw peak labels
label_fun <- function(x) {
    ints <- unlist(intensity(x))
    mzs <- format(unlist(mz(x)), digits = 4)
    mzs[ints < 5] <- ""
    mzs
}

##-----------------------------------------------------------------
## Pre-process MS2 spectra
##-----------------------------------------------------------------

#' processing on spectra with one precursor mass
# inputs: 
# x is precursor mass, 
# spec is the spectra file (sps_all is mzML input processed spectra, gnps, hmdb or mbank), 
# ppmx is ppm value
spec2_Processing <- function(z, obj, spec = "spec_all", ppmx = 15){
    if (spec == "spec_all"){
        #' Subset the dataset to MS2 spectra matching the m/z
        sps <- filterPrecursorMzRange(obj, mz = z + ppm(c(-z, z), 10))
    } else if (spec == "gnps"){
        #gnps spectra that contains precursor mass
        has_mz <- containsMz(obj, mz = z, ppm = ppmx)
        #' Subset the GNPS Spectra
        sps <- obj[has_mz]
    } else if (spec == "hmdb"){
        #hmdb spectra that contains precursor mass
        has_mz <- containsMz(obj, mz = z, ppm = ppmx)
        #' Subset the HMDB Spectra
        sps <- obj[has_mz]
    } else if (spec == "mbank"){
        has_mz <- containsMz(obj, mz = z, ppm = ppmx)
        #' Subset the HMDB Spectra
        sps <- obj[has_mz]
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
        removePrecursorPeaks <- function(m){
            m <- m[m[, "mz"] < z, ]
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
# Usage:
# spec2_Processing(x = 231.15, spec = "sps_all", ppmx = 15)


##-----------------------------------------------------------------
## Extract peaksdata in a dataframe
##-----------------------------------------------------------------

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

# ---------- spec_dereplication ----------

spec_dereplication<- function(pre_tbl, proc_mzml, db, result_dir, file_id, input_dir, ppmx, error = TRUE){
    
    
    ####-------------------------------------------------------------
    #### Dereplication with all or GNPS ----
    ####-------------------------------------------------------------

    
    sps_all <- Spectra(proc_mzml, backend = MsBackendMzR())
        
    tbl <- read.table(pre_tbl)
    pre_mz <- tbl[[1]]
    
    if (db == "all" || db =="gnps"){

        load(file = paste(input_dir,"gnps.rda", sep = ""))
        
        print(R.version)
        print(packageVersion("Spectra"))
        load(file = paste(input_dir,"gnps.rda", sep = ""))
        
        # common

        id_X <- c()
        premz <- c()
        rtmin <- c()
        rtmax <- c()
        rtmed <- c()
        rtmean <- c()

        # gnps
        GNPSmax_similarity <- c()
        GNPSmzScore <- c()
        GNPSintScore <- c()
        GQMatchingPeaks <- c()
        GNPSTotalPeaks <- c()
        gQueryTotalPeaks<- c()
        GNPSSMILES <- c()
        GNPSspectrumID <- c()
        GNPScompound_name <- c()
        GNPSmirrorSpec <- c()
        
        # common
        Source <- c()
        
        nx <- 0
        

        for (x in pre_mz){

            
            nx <- nx+1
            
            spsrt <- filterPrecursorMzRange(sps_all, x)
        
            
            id_Xx <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                            "R", as.character(round(median(spsrt$rtime, na.rm = TRUE), digits = 0)), 
                            "ID", as.character(nx), sep = '')
            id_X <- c(id_X, id_Xx)

            pre <- x
            premz <- c(premz, pre)

            rti <- min(spsrt$rtime)
            rtmin <- c(rtmin, rti)

            rtx <- max(spsrt$rtime)
            rtmax <- c(rtmax, rtx)


            rtmd <- median(spsrt$rtime, na.rm = TRUE)
            rtmed <- c(rtmed, rtmd)

            rtmn <- mean(spsrt$rtime, na.rm = TRUE)
            rtmean <- c(rtmean, rtmn)
       
            #### input spec with pre_mz
            sps <- spec2_Processing(x, sps_all, spec = "spec_all")
        
            #### GNPS spec with pre_mz
            gnps_with_mz <- spec2_Processing(x, gnpsdb, spec = "gnps", ppmx) # change here later

        
            dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/GNPS/", sep = ""), "./"), sep ="")
            if (!file.exists(dir_name)){
                dir.create(dir_name, recursive = TRUE)
            }
        
            if (length(sps) > 1 && length(gnps_with_mz) >1){
                #' Compare experimental spectra against GNPS
                res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
            
                #' obtain GNPS spectra that matches the most with m/z MS2 spectra
                idx <- which(res == max(res), arr.ind = TRUE)
                
                if (nrow(idx)>1){
                    idx <- idx[1,]
                }
        
                
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
                
                    GNPSscore <- max(res)
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)
                
                    GNPSmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)
            
                
                    GNPSint <- mean(1-(df_peaklists[,"diff"]/100))
                    GNPSintScore <- c(GNPSintScore, GNPSint)
                    
                
                    GQMatPeaks <- nrow(df_peaklists)
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)
                    
                    GNPSTPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)
                
                    gQTPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)
                
                
                    GNPS_SMILES <- gnps_best_match$SMILES
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)
                
                
                    GNPSID <- gnps_best_match$SPECTRUMID
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)
                
                    GNPSname <- gnps_best_match$NAME
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)
                
                
                    GNPSSpec <- str_replace(name_plotmirror, input_dir, "./")
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)
                
                    Src <- "GNPS"
                    Source <- c(Source, Src)
                }
                else{

            
                    GNPSscore <- NA
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)
                    
                    GNPSmz <- NA
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)
                

                    GNPSint <- NA
                    GNPSintScore <- c(GNPSintScore, GNPSint)
                    
                
                    GQMatPeaks <- NA
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)
                
                
                    GNPSTPeaks <- NA
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)
                
                
                    gQTPeaks<- NA
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)
                
                
                    GNPS_SMILES <- NA
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)
                
                
                    GNPSID <- NA
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)
                
                
                    GNPSname <- NA
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)
                
                
                    GNPSSpec <- NA
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)
                
                
                    Src <- NA
                    Source <- c(Source, Src)
            
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

            
                    GNPSscore <- max(res)
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)
                    
            
                    GNPSmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)))
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)
                
                
                    GNPSint <- mean(1-(df_peaklists[,"diff"]/100))
                    GNPSintScore <- c(GNPSintScore, GNPSint)
                
                
                    GQMatPeaks <- nrow(df_peaklists)
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)
                

                    GNPSTPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)
                
             
                    gQTPeaks<- nrow(peaksData(sps))
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)
                
                    GNPS_SMILES <- gnps_best_match$SMILES
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)
                    
                    GNPSID <- gnps_best_match$SPECTRUMID
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)
                    
                    GNPSname <- gnps_best_match$NAME
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)
            
                
                    GNPSSpec <- str_replace(name_plotmirror, input_dir, "./")
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)
                
                
                    Src <- "GNPS"
                    Source <- c(Source, Src)
                }
                else{

            
                    GNPSscore <- NA
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)
                
                    GNPSmz <- NA
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)
                
                
                    GNPSint <- NA
                    GNPSintScore <- c(GNPSintScore, GNPSint)
                
                
                    GQMatPeaks <- NA
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)
                
                
                    GNPSTPeaks <- NA
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)
                
                
                    gQTPeaks<- NA
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)
                
                
                    GNPS_SMILES <- NA
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)
                
                
                    GNPSID <- NA
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)
                
                
                    GNPSname <- NA
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)
                
                
                    GNPSSpec <- NA
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)
                
                
                    Src <- NA
                    Source <- c(Source, Src)
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
                
                    GNPSscore <- max(res)
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)
                    
                    GNPSmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)))
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)
                
                    GNPSint <- mean(1-(df_peaklists[,"diff"]/100))
                    GNPSintScore <- c(GNPSintScore, GNPSint)
                
                
                
                    GQMatPeaks <- nrow(df_peaklists)
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)
                
                
                    GNPSTPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)
                

                    gQTPeaks<- nrow(peaksData(sps))
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)
                
                    GNPS_SMILES <- gnps_best_match$SMILES
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)
                
                    GNPSID <- gnps_best_match$SPECTRUMID
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)
                
                
                    GNPSname <- gnps_best_match$NAME
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)
                
                    
                    GNPSSpec <- str_replace(name_plotmirror, input_dir, "./")
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)
                
                
                    Src <- "GNPS"
                    Source <- c(Source, Src)
                    }
            
                else{

                    GNPSscore <- NA
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)
                
                
                    GNPSmz <- NA
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)
                
                
                    GNPSint <- NA
                    GNPSintScore <- c(GNPSintScore, GNPSint)
                
                
                    GQMatPeaks <- NA
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)
                
                
                    GNPSTPeaks <- NA
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)
                
                
                    gQTPeaks<- NA
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)
                    
                
                    GNPS_SMILES <- NA
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)
                    
                
                    GNPSID <- NA
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)
                
                
                    GNPSname <- NA
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)
                
                
                    GNPSSpec <- NA
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)
                
                
                    Src <- NA
                    Source <- c(Source, Src)
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

            
                    GNPSscore <- max(res)
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)
                
                
                
                    GNPSmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)
                
                    GNPSint <- mean(1-(df_peaklists[,"diff"]/100))
                    GNPSintScore <- c(GNPSintScore, GNPSint)
                
                    GQMatPeaks <- NA
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)
                
            
                    GNPSTPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)                
                
                    gQTPeaks<- nrow(peaksData(sps)[[1]])
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)


                    GNPS_SMILES <- gnps_best_match$SMILES
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)

                    GNPSname <- gnps_best_match$NAME
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)


                    GNPSID <- gnps_best_match$SPECTRUMID
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)


                    GNPSSpec <- NA
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)


                    Src <- "GNPS"
                    Source <- c(Source, Src)
                }
                else{


                    GNPSscore <- NA
                    GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)


                    GNPSmz <- NA
                    GNPSmzScore <- c(GNPSmzScore, GNPSmz)


                    GNPSint <- NA
                    GNPSintScore <- c(GNPSintScore, GNPSint)


                    GQMatPeaks <- NA
                    GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)


                    GNPSTPeaks <- NA
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)


                    gQTPeaks<- NA
                    gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)


                    GNPS_SMILES <- NA
                    GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)


                    GNPSID <- NA
                    GNPSspectrumID <- c(GNPSspectrumID, GNPSID)


                    GNPSname <- NA
                    GNPScompound_name <- c(GNPScompound_name, GNPSname)


                    GNPSSpec <- NA
                    GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)


                    Src <- NA
                    Source <- c(Source, Src)
                }

            }
            else{


                GNPSscore <- NA
                GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)

                GNPSmz <- NA
                GNPSmzScore <- c(GNPSmzScore, GNPSmz)


                GNPSint <- NA
                GNPSintScore <- c(GNPSintScore, GNPSint)


                GQMatPeaks <- NA
                GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)


                GNPSTPeaks <- NA
                GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)


                gQTPeaks<- NA
                gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)


                GNPS_SMILES <- NA
                GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)


                GNPSID <- NA
                GNPSspectrumID <- c(GNPSspectrumID, GNPSID)

                GNPSname <- NA
                GNPScompound_name <- c(GNPScompound_name, GNPSname)


                GNPSSpec <- NA
                GNPSmirrorSpec <- c(GNPSmirrorSpec, GNPSSpec)


                Src <- NA
                Source <- c(Source, Src)

            }
            
        }
        df_gnps <- cbind(id_X, premz, rtmin, rtmax, rtmed, rtmean, GNPSmax_similarity, GNPSmzScore, 
                        GNPSintScore, GQMatchingPeaks, GNPSTotalPeaks, gQueryTotalPeaks, 
                        GNPSSMILES, GNPSspectrumID, GNPScompound_name, GNPSmirrorSpec, Source)

        #write.csv(df_gnps, str_remove(paste(input_dir, result_dir, "/spectral_dereplication/gnps.csv", sep = ""), "./"))
        write.csv(df_gnps, paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/gnps.csv", sep = ""), "./"), sep = ""))
    }
    
    ####-------------------------------------------------------------
    #### Dereplication with all or HMDB ----
    ####-------------------------------------------------------------
    if (db == "all" || db =="hmdb"){
<<<<<<< HEAD
        
        load(file = paste(input_dir,"hmdb.rda", sep = ""))
        
=======
        load(file = paste(input_dir,"hmdb.rda", sep = ""))
>>>>>>> 25c6491 (cleaned directory)
        # common

        id_X <- c()
        premz <- c()
        rtmin <- c()
        rtmax <- c()
        rtmed <- c()
        rtmean <- c()

        # hmdb 
        HMDBmax_similarity <- c()
        HMDBmzScore <- c()
        HMDBintScore <- c()
        HQMatchingPeaks <- c()
        HMDBTotalPeaks <- c()
        hQueryTotalPeaks<- c()
        HMDBcompoundID <- c()
        HMDBmirrorSpec <- c()
    
        # common
        Source <- c()
        nx <- 0
        for (x in pre_mz){
            
            nx <- nx+1
            
            spsrt <- filterPrecursorMzRange(sps_all, x)
        
            
            id_Xx <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                            "R", as.character(round(median(spsrt$rtime, na.rm = TRUE), digits = 0)), 
                            "ID", as.character(nx), sep = '')
            id_X <- c(id_X, id_Xx)

            pre <- x
            premz <- c(premz, pre)

            rti <- min(spsrt$rtime)
            rtmin <- c(rtmin, rti)

            rtx <- max(spsrt$rtime)
            rtmax <- c(rtmax, rtx)


            rtmd <- median(spsrt$rtime, na.rm = TRUE)
            rtmed <- c(rtmed, rtmd)

            rtmn <- mean(spsrt$rtime, na.rm = TRUE)
            rtmean <- c(rtmean, rtmn)

            #### input spec with pre_mz
            sps <- spec2_Processing(x, sps_all, spec = "spec_all", ppmx = NULL)

            #### HMDB spec with pre_mz
            hmdb_with_mz <- spec2_Processing(x, hmdb, spec = "hmdb", ppmx = 15)


            dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/HMDB/", sep = ""), "./"), sep ="")
            if (!file.exists(dir_name)){
                dir.create(dir_name, recursive = TRUE)
            }
            if (length(sps) > 1 && length(hmdb_with_mz) > 1){
                #' Compare experimental spectra against HMDB
                res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)

                #' obtain HMDB spectra that matches the most with m/z MS2 spectra
                idx <- which(res == max(res), arr.ind = TRUE)
                
                if (nrow(idx)>1){
                    idx <- idx[1,]
                }
                
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


                    HMDBscore <- max(res)
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)

                    HMDBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)

                    HMDBint <- mean(1-(df_peaklists[,"diff"]/100))
                    HMDBintScore <- c(HMDBintScore, HMDBint)


                    HQMatPeaks <- nrow(df_peaklists)
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)

                    hQTPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- hmdb_best_match$compound_id
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)

                    HMDBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)

                    Src <- "HMDB"
                    Source <- c(Source, Src)
                }
                else{

                    HMDBscore <- NA
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                    HMDBmz <- NA
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)


                    HMDBint <- NA
                    HMDBintScore <- c(HMDBintScore, HMDBint)


                    HQMatPeaks <- NA
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- NA
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                    hQTPeaks<- NA
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- "none"
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)


                    HMDBSpec <- NA
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                    Src <- NA
                    Source <- c(Source, Src)
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

                    HMDBscore <- max(res)
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                    HMDBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)

                    HMDBint <- mean(1-(df_peaklists[,"diff"]/100))
                    HMDBintScore <- c(HMDBintScore, HMDBint)


                    HQMatPeaks <- nrow(df_peaklists)
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                    hQTPeaks<- nrow(peaksData(sps)[[1]])
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- hmdb_best_match$compound_id
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)



                    HMDBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                    Src <- "HMDB"
                    Source <- c(Source, Src)
                }
                else{

                    HMDBscore <- NA
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                    HMDBmz <- NA
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)


                    HMDBint <- NA
                    HMDBintScore <- c(HMDBintScore, HMDBint)


                    HQMatPeaks <- NA
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- NA
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                    hQTPeaks<- NA
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- "none"
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)


                    HMDBSpec <- NA
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                    Src <- NA
                    Source <- c(Source, Src)
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

                    HMDBscore <- max(res)
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                    HMDBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)

                    HMDBint <- mean(1-(df_peaklists[,"diff"]/100))
                    HMDBintScore <- c(HMDBintScore, HMDBint)


                    HQMatPeaks <- nrow(df_peaklists)
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                    hQTPeaks<- nrow(peaksData(sps)[[1]])
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- hmdb_best_match$compound_id
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)



                    HMDBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                    Src <- "HMDB"
                    Source <- c(Source, Src)
                    }
                else{

                    HMDBscore <- NA
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                    HMDBmz <- NA
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)


                    HMDBint <- NA
                    HMDBintScore <- c(HMDBintScore, HMDBint)



                    HQMatPeaks <- NA
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- NA
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                    hQTPeaks<- NA
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- "none"
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)


                    HMDBSpec <- NA
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                    Src <- NA
                    Source <- c(Source, Src)
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

                    HMDBscore <- max(res)
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                    HMDBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)


                    HMDBint <- mean(1-(df_peaklists[,"diff"]/100))
                    HMDBintScore <- c(HMDBintScore, HMDBint)



                    HQMatPeaks <- nrow(df_peaklists)
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                    hQTPeaks<- nrow(peaksData(sps)[[1]])
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- hmdb_best_match$compound_id
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)



                    HMDBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                    Src <- "HMDB"
                    Source <- c(Source, Src)
                }
                else{

                    HMDBscore <- NA
                    HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                    HMDBmz <- NA
                    HMDBmzScore <- c(HMDBmzScore, HMDBmz)


                    HMDBint <- NA
                    HMDBintScore <- c(HMDBintScore, HMDBint)


                    HQMatPeaks <- NA
                    HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                    HMDBTPeaks <- NA
                    HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                    hQTPeaks<- NA
                    hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                    HMDBID <- "none"
                    HMDBcompoundID <- c(HMDBcompoundID, HMDBID)


                    HMDBSpec <- NA
                    HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                    Src <- NA
                    Source <- c(Source, Src)
                }

            }
            else{


                HMDBscore <- NA
                HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)


                HMDBmz <- NA
                HMDBmzScore <- c(HMDBmzScore, HMDBmz)


                HMDBint <- NA
                HMDBintScore <- c(HMDBintScore, HMDBint)


                HQMatPeaks <- NA
                HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                HMDBTPeaks <- NA
                HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)


                hQTPeaks<- NA
                hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                HMDBID <- "none"
                HMDBcompoundID <- c(HMDBcompoundID, HMDBID)

                HMDBSpec <- NA
                HMDBmirrorSpec <- c(HMDBmirrorSpec, HMDBSpec)


                Src <- NA
                Source <- c(Source, Src)
            }
        }
        df_hmdb <- cbind(id_X, premz, rtmin, rtmax, rtmed, rtmean, HMDBmax_similarity, HMDBmzScore, 
                         HMDBintScore, HQMatchingPeaks, HMDBTotalPeaks, hQueryTotalPeaks, 
                         HMDBcompoundID, HMDBmirrorSpec, Source)

        #write.csv(df_hmdb, str_remove(paste(result_dir, "/spectral_dereplication/hmdb.csv", sep = ""), "./"))
        write.csv(df_hmdb, paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/hmdb.csv", sep = ""), "./"), sep = ""))
    }
    ####-------------------------------------------------------------
    #### Dereplication with all or MassBank ----
     ####-------------------------------------------------------------
    if (db == "all" || db =="mbank"){
<<<<<<< HEAD

        load(file = paste(input_dir,"mbankNIST.rda", sep = ""))
=======
        
        load(file = paste(input_dir,"mbank.rda", sep = ""))
>>>>>>> 25c6491 (cleaned directory)
        
        # common

        id_X <- c()
        premz <- c()
        rtmin <- c()
        rtmax <- c()
        rtmed <- c()
        rtmean <- c()

        # mbank
        MBmax_similarity <- c()
        MBmzScore <- c()
        MBintScore <- c()
        MQMatchingPeaks <- c()
        MBTotalPeaks <- c()
        mQueryTotalPeaks<- c()
        MBformula <- c()
        MBinchiKEY <- c()
        MBspectrumID <- c()
        MBcompound_name <- c()
        MBmirrorSpec <- c()

        # common
        Source <- c()
        nx <- 0
        for (x in pre_mz){

            
            nx <- nx+1
            
            
            spsrt <- filterPrecursorMzRange(sps_all, x)
        
            
            id_Xx <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                            "R", as.character(round(median(spsrt$rtime, na.rm = TRUE), digits = 0)), 
                            "ID", as.character(nx), sep = '')
            id_X <- c(id_X, id_Xx)

            pre <- x
            premz <- c(premz, pre)

            rti <- min(spsrt$rtime)
            rtmin <- c(rtmin, rti)

            rtx <- max(spsrt$rtime)
            rtmax <- c(rtmax, rtx)


            rtmd <- median(spsrt$rtime, na.rm = TRUE)
            rtmed <- c(rtmed, rtmd)

            rtmn <- mean(spsrt$rtime, na.rm = TRUE)
            rtmean <- c(rtmean, rtmn)

            #### input spec with pre_mz
            sps <- spec2_Processing(x, sps_all, spec = "spec_all")

            #### GNPS spec with pre_mz
            mbank_with_mz <- spec2_Processing(x, mbank, spec = "mbank", ppmx = 15)
            
            

            dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/MassBank/", sep = ""), "./"), sep ="")
            if (!file.exists(dir_name)){
                dir.create(dir_name, recursive = TRUE)
            }
            if (length(sps) > 1 && length(mbank_with_mz) >1){

                #' Compare experimental spectra against MassBank
                res <- compareSpectra(sps, mbank_with_mz, ppm = 15)

                #' obtain GNPS spectra that matches the most with m/z MS2 spectra
                idx <- which(res == max(res), arr.ind = TRUE)
                
                if (nrow(idx)>1){
                    idx <- idx[1,]
                }
                
                
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

                    mbscore <- max(res)
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- mean(1-(df_peaklists[,"diff"]/100))
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- nrow(df_peaklists)
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks) 

                    MBTPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- nrow(peaksData(sps[idx[1]])[[1]]) 
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)


                    MBfor <- mbank_best_match$formula
                    MBformula<- c(MBformula, MBfor)

                    MBinchiK <- mbank_best_match$inchikey
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- mbank_best_match$accession
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- mbank_best_match$name
                    MBcompound_name <- c(MBcompound_name, MBname)


                    MBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- "MassBank"
                    Source <- c(Source, Src)
                }
                else{

                    mbscore <- NA
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- NA
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- NA
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- NA
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                    MBTPeaks <- NA
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- NA 
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                    MBfor <- NA
                    MBformula<- c(MBformula, MBfor)

                    MBinchiK <- NA
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- NA
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- NA
                    MBcompound_name <- c(MBcompound_name, MBname)

                    MBSpec <- NA
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- NA
                    Source <- c(Source, Src)
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

                    mbscore <- max(res)
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- mean(1-(df_peaklists[,"diff"]/100))
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- nrow(df_peaklists)
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                    MBTPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- nrow(peaksData(sps)[[1]])
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                    MBfor <- mbank_best_match$formula
                    MBformula<- c(MBformula, MBfor)

                    MBinchiK <- mbank_best_match$inchikey
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- mbank_best_match$accession
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- mbank_best_match$name
                    MBcompound_name <- c(MBcompound_name, MBname)

                    MBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- "MassBank"
                    Source <- c(Source, Src)
                }
                else{

                    mbscore <- NA
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- NA
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- NA
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- NA
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                    MBTPeaks <- NA
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- NA
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                    MBfor <- NA
                    MBformula<- c(MBformula, MBfor)

                    MBinchiK <- NA
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- NA
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- NA
                    MBcompound_name <- c(MBcompound_name, MBname)

                    MBSpec <- NA
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- NA
                    Source <- c(Source, Src)
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

                    mbscore <- max(res)
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- mean(1-(df_peaklists[,"diff"]/100))
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- nrow(df_peaklists)
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)


                    MBTPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- nrow(peaksData(sps)[[1]])
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                    MBfor <- mbank_best_match$formula
                    MBformula<- c(MBformula, MBfor)

                    MBinchiK <- mbank_best_match$inchikey
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- mbank_best_match$accession
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- mbank_best_match$name
                    MBcompound_name <- c(MBcompound_name, MBname)

                    MBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- "MassBank"
                    Source <- c(Source, Src)
                    }
                else{

                    mbscore <- NA
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- NA
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- NA
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- NA
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                    MBTPeaks <- NA
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- NA
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                    MBfor <- NA
                    MBformula<- c(MBformula, MBfor)


                    MBinchiK <- NA
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- NA
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- NA
                    MBcompound_name <- c(MBcompound_name, MBname)

                    MBSpec <- NA
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- NA
                    Source <- c(Source, Src)
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


                    mbscore <- max(res)
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- mean(1-(df_peaklists[,"diff"]/100))
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- nrow(df_peaklists)
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                    MBTPeaks <- nrow(peaksData(mbank_best_match)[[1]])
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- nrow(peaksData(sps)[[1]])
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                    MBfor <- mbank_best_match$formula
                    MBformula<- c(MBformula, MBfor)

                    MBinchiK <- mbank_best_match$inchikey
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- mbank_best_match$accession
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- mbank_best_match$name
                    MBcompound_name <- c(MBcompound_name, MBname)

                    MBSpec <- str_replace(name_plotmirror, input_dir, "./")
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- "MassBank"
                    Source <- c(Source, Src)
                }
                else{

                    mbscore <- NA
                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                    MBmz <- NA
                    MBmzScore <- c(MBmzScore, MBmz)

                    MBint <- NA
                    MBintScore <- c(MBintScore, MBint)

                    MQMatPeaks <- NA
                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                    MBTPeaks <- NA
                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                    mQTPeaks<- NA
                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                    MBfor <- NA
                    MBformula<- c(MBformula, MBfor)


                    MBinchiK <- NA
                    MBinchiKEY <- c(MBinchiKEY, MBinchiK)

                    MBID <- NA
                    MBspectrumID<- c(MBspectrumID, MBID)

                    MBname <- NA
                    MBcompound_name <- c(MBcompound_name, MBname)

                    MBSpec <- NA
                    MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                    Src <- NA
                    Source <- c(Source, Src)
                }

            }
            else{

                mbscore <- NA
                MBmax_similarity<- c(MBmax_similarity, mbscore)

                MBmz <- NA
                MBmzScore <- c(MBmzScore, MBmz)

                MBint <- NA
                MBintScore <- c(MBintScore, MBint)

                MQMatPeaks <- NA
                MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                MBTPeaks <- NA
                MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                mQTPeaks<- NA
                mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)

                MBfor <- NA
                MBformula<- c(MBformula, MBfor)

                MBinchiK <- NA
                MBinchiKEY <- c(MBinchiKEY, MBinchiK)


                MBID <- NA
                MBspectrumID<- c(MBspectrumID, MBID)

                MBname <- NA
                MBcompound_name <- c(MBcompound_name, MBname)

                MBSpec <- NA
                MBmirrorSpec <- c(MBmirrorSpec, MBSpec)

                Src <- NA
                Source <- c(Source, Src)

            }
        }
        
        df_mbank <- cbind(id_X, premz, rtmin, rtmax, rtmed, rtmean, MBmax_similarity, MBmzScore, MBintScore, MQMatchingPeaks, 
                         MBTotalPeaks, mQueryTotalPeaks, MBformula, MBinchiKEY, MBspectrumID, MBcompound_name, MBmirrorSpec, 
                         Source)
        write.csv(df_mbank, paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/mbank.csv", sep = ""), "./"), sep = ""))


    }
    
}
spec_dereplication(pre_tbl, proc_mzml, db, result_dir, file_id, input_dir, ppmx)

