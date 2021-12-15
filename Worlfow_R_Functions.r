##-----------------------------------------------------------------
## filter intensity 
##-----------------------------------------------------------------

#' Define a filtering function and remove peaks less than 0.05 of intensity
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * 0.05
}

# Usage:
# filterIntensity(spectra_object, intensity = low_int)

#' filterIntensity is a predefined function in Spectra package

##-----------------------------------------------------------------
## normalize intensity 
##-----------------------------------------------------------------

#' Define a function to *normalize* the intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}

# Usage:
# addProcessing(sps, norm_int)

#' addProcessing is a predefined function in Spectra package

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
# Usage: 
# plotSpectraMirror(sps[idx[1]], gnps_with_mz[idx[2]], tolerance = 0.2,
# labels = label_fun, labelPos = 2, labelOffset = 0.2, labelSrt = -30)

#' plotSpectraMirror is a predefined function in Spectra package

##-----------------------------------------------------------------
## Result Directories 
##-----------------------------------------------------------------

#' Specifying a function for creating result directories for each input mzml
# inputs for the function:
# x is mzML File, input directory and working directory
ms2_rfilename<- function(input_dir){
    #list_ms2_files <- intersect(list.files(input_dir, pattern = "_PRM_"), list.files(input_dir, pattern = ".mzML"))
    list_ms2_files <- list.files(input_dir, pattern = ".mzML")
    mzml_files <- paste(input_dir, list_ms2_files, sep = "")
    
    #store the result file names to return to this function as output
    ResultFileNames <- c()
    File_id <- c()
    nx <- 0
    # x is mzML files
    for (i in 1:length(mzml_files)){
        nx <- nx+1
        # remove .mzML to extract just the names
        name_mzml <- str_remove(as.character(mzml_files[i]), ".mzML")      
        #' for each file a subdirectory is created to store all results in that, add working directory
        if (!file.exists(name_mzml)){
            dir.create(name_mzml) ##create folder
        }
        ResultFileNames<- c(ResultFileNames, name_mzml)
        File_id <- c(File_id, paste("file_", nx, sep = ""))
    }
    input_table <- cbind(mzml_files, ResultFileNames, File_id)
    return(input_table)
}

#' All spectra in mzML files preprocessing, return two outputs, pre-processed MS2 spectra and all precursor masses
# x is one mzML file
spec_Processing <- function(x){
    # read the spectra
    sps_all <- Spectra(x, backend = MsBackendMzR())
    if (length(na.omit(unique(sps_all$precursorMz)))>1){
        #' Change backend to a MsBackendDataFrame: load data into memory
        sps_all <- setBackend(sps_all, MsBackendDataFrame())
        #' Filter Empty Spectra
        sps_all <- filterEmptySpectra(sps_all)
        #' Extract Precursor m/z(s) in each file
        pre_mz <- unique(precursorMz(sps_all))
        #' Remove any NAs
        pre_mz <- na.omit(pre_mz)
        spsall_pmz <- list(sps_all, pre_mz)
        return(spsall_pmz)
    }
    else{
        #' Filter Empty Spectra
        sps_all <- filterEmptySpectra(sps_all)
        #' Extract Precursor m/z(s) in each file
        pre_mz <- unique(precursorMz(sps_all))
        #' Remove any NAs
        pre_mz <- na.omit(pre_mz)
        spsall_pmz <- list(sps_all, pre_mz)
        return(spsall_pmz)
    }
}

#' processing on spectra with one precursor mass
# inputs: x is precursor mass, spec is the spectra file (sps_all is mzML input processed spectra, gnps, hmdb), ppmx is ppm value
spec2_Processing <- function(x, spec = "sps_all", ppmx = 15){
    if (spec == "sps_all"){
        #' Subset the dataset to MS2 spectra matching the m/z
        sps <- filterPrecursorMz(sps_all, mz = x + ppm(c(-x, x), 10))
    } else if (spec == "gnps"){
        #gnps spectra that contains precursor mass
        has_mz <- containsMz(gnps, mz = x, ppm = ppmx)
        #' Subset the GNPS Spectra
        sps <- gnps[has_mz]
    } else if (spec == "hmdb"){
        #hmdb spectra that contains precursor mass
        has_mz <- containsMz(hmdb, mz = x, ppm = ppmx)
        #' Subset the HMDB Spectra
        sps <- hmdb[has_mz]
    } else if (spec == "mbank"){
        has_mz <- containsMz(mbank,mz = x, ppm = ppmx)
        #' Subset the HMDB Spectra
        sps <- mbank[has_mz]
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
            y <- y[y[, "mz"] < x, ]
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
        # if 0 peaks, remove the relevant spectra from gnps or hmdb
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
peakdf <- function(a, b){
    ##----------------------------------------------------------------
    ## SCORING -
    ##----------------------------------------------------------------
    #' obtain peaklists for both query spectra and best matched spectra from MassBank
    z<- peaksData(a)[[1]] #from GNPS/HMDB
    y <- peaksData(b)[[1]] #from query
    if (!(nrow(z)==0)){
    #' Since we used 15 ppm, so to find the range, calculate the mass error range
    range <-  (15/1000000)*y[,"mz"]
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


# x is one pre_mz, db is GNPS, HMDB, MassBank
spec_dereplication<- function(x, db, result_dir, file_id, input_dir){
    if (db == "GNPS"){
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(a, spec = "sps_all", ppmx = NULL)
        
        #### GNPS spec with pre_mz
        gnps_with_mz <- spec2_Processing(a, spec = "gnps", ppmx = 15)
        
        
        if (length(sps) > 1 && length(gnps_with_mz) >1){
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            gnps_best_match <- gnps_with_mz[idx[2]]
            df_peaklists <- peakdf(gnps_best_match, sps[idx[1]])
            
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/GNPS/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            gnps_best_match <- gnps_with_mz[gx]

            df_peaklists <- peakdf(gnps_best_match, sps)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/GNPS/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
            #' obtain MB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(gnps_with_mz, sps[gx])
            gnps_best_match <- gnps_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/GNPS/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
            gnps_best_match <- gnps_with_mz
            df_peaklists <- peakdf(gnps_best_match, sps)
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/GNPS/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
        
        df_gnps <- cbind(id_X, premz, rtmin, rtmax, rtmed, GNPSmax_similarity, GNPSmzScore, 
                        GNPSintScore, GQMatchingPeaks, GNPSTotalPeaks, gQueryTotalPeaks, 
                        GNPSSMILES, GNPSspectrumID, GNPScompound_name, GNPSmirrorSpec, Source)
        return(df_gnps)
    }
    if (db == "HMDB"){
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(a, spec = "sps_all", ppmx = NULL)
        
        #### GNPS spec with pre_mz
        hmdb_with_mz <- spec2_Processing(a, spec = "hmdb", ppmx = 15)
        
        if (length(sps) > 1 && length(hmdb_with_mz) >1){
            #' Compare experimental spectra against HMDB
            res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            hmdb_best_match <- hmdb_with_mz[idx[2]]
            df_peaklists <- peakdf(hmdb_best_match, sps[idx[1]])
            
            if (!(is.null(df_peaklists))){
                print("more sps and more hmsb_with_mz")
                dir_name <- paste(result_dir, "/spectral_dereplication/HMDB/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
                 
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                HMDBcompoundID <- gnps_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
            }
            else{
                print("NO more sps and more hmsb_with_mz")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- NA
                HMDBmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(hmdb_with_mz) >1){
            
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            hmdb_best_match <- hmdb_with_mz[gx]

            df_peaklists <- peakdf(hmdb_best_match, sps)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                print("one sps and more hmsb_with_mz")
                dir_name <- paste(result_dir, "/spectral_dereplication/HMDB/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
                 
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                HMDBcompoundID <- gnps_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
            }
            else{
                print("NO one sps and more hmsb_with_mz")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- NA
                HMDBmirrorSpec <- NA
                Source <- NA
            }
        
        }
        else if (length(sps) > 1 && length(hmdb_with_mz) == 1){
            #' obtain MB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(hmdb_with_mz, sps[gx])
            hmdb_best_match <- hmdb_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                print("more sps and one hmsb_with_mz")
                dir_name <- paste(result_dir, "/spectral_dereplication/HMDB/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
                 
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                HMDBcompoundID <- gnps_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
                }
            else{
                print("NO more sps and one hmsb_with_mz")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- NA
                HMDBmirrorSpec <- NA
                Source <- NA
            }
        }
        else if (length(sps) == 1 && length(hmdb_with_mz) == 1){
            hmdb_best_match <- hmdb_with_mz
            df_peaklists <- peakdf(hmdb_best_match, sps)
            if (!(is.null(df_peaklists))){
                print("one sps and one hmsb_with_mz")
                dir_name <- paste(result_dir, "/spectral_dereplication/HMDB/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
                 
                HMDBmax_similarity <- max(res)
                HMDBmzScore <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                HMDBintScore <- mean(1-(df_peaklists[,"diff"]/100))
                HQMatchingPeaks <- nrow(df_peaklists)
                HMDBTotalPeaks <- nrow(peaksData(hmdb_best_match)[[1]])
                hQueryTotalPeaks<- nrow(peaksData(sps[idx[1]])[[1]])
                HMDBcompoundID <- gnps_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
            }
            else{
                print("NO one sps and one hmsb_with_mz")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                HMDBmax_similarity <- NA
                HMDBmzScore <- NA
                HMDBintScore <- NA
                HQMatchingPeaks <- NA
                HMDBTotalPeaks <- NA
                hQueryTotalPeaks<- NA
                HMDBcompoundID <- NA
                HMDBmirrorSpec <- NA
                Source <- NA
            }
                
        }
        else{
            print("NO sps and NO hmsb_with_mz")
            id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
            premz <- x
            rtmin <- max(sps$rtime)
            rtmax <- min(sps$rtime)
            rtmed <- median(sps$rtime, na.rm = TRUE)
            HMDBmax_similarity <- NA
            HMDBmzScore <- NA
            HMDBintScore <- NA
            HQMatchingPeaks <- NA
            HMDBTotalPeaks <- NA
            hQueryTotalPeaks<- NA
            HMDBcompoundID <- NA
            HMDBmirrorSpec <- NA
            Source <- NA
        }
        
        df_hmdb <- cbind(id_X, premz, rtmin, rtmax, rtmed, HMDBmax_similarity, HMDBmzScore, 
                         HMDBintScore, HQMatchingPeaks, HMDBTotalPeaks, hQueryTotalPeaks, 
                         HMDBcompoundID, HMDBmirrorSpec, Source)
        return(df_hmdb)
    }
    if (db == "MB"){
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(a, spec = "sps_all", ppmx = NULL)
        
        #### GNPS spec with pre_mz
        mbank_with_mz <- spec2_Processing(a, spec = "mbank", ppmx = 15)
        
        
        if (length(sps) > 1 && length(mbank_with_mz) >1){
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, mbank_with_mz, ppm = 15)
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            mbank_best_match <- mbank_with_mz[idx[2]]
            df_peaklists <- peakdf(mbank_best_match, sps[idx[1]])
            
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/MB/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            gnps_best_match <- gnps_with_mz[gx]

            df_peaklists <- peakdf(gnps_best_match, sps)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/GNPS/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
            #' obtain MB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(gnps_with_mz, sps[gx])
            gnps_best_match <- gnps_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/GNPS/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
            gnps_best_match <- gnps_with_mz
            df_peaklists <- peakdf(gnps_best_match, sps)
            if (!(is.null(df_peaklists))){
                dir_name <- paste(result_dir, "/spectral_dereplication/GNPS/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
        
        df_gnps <- cbind(id_X, premz, rtmin, rtmax, rtmed, GNPSmax_similarity, GNPSmzScore, 
                        GNPSintScore, GQMatchingPeaks, GNPSTotalPeaks, gQueryTotalPeaks, 
                        GNPSSMILES, GNPSspectrumID, GNPScompound_name, GNPSmirrorSpec, Source)
        return(df_gnps)
    }
    
}

save.image(file = "R_Functions.RData")

































save.image(file = "R_Functions.RData")




