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
    #' Change backend to a MsBackendDataFrame: load data into memory
    #sps_all <- setBackend(sps_all, MsBackendDataFrame())
    #' Filter Empty Spectra
    sps_all <- filterEmptySpectra(sps_all)
    #' Extract Precursor m/z(s) in each file
    pre_mz <- unique(precursorMz(sps_all))
    #' Remove any NAs
    pre_mz <- na.omit(pre_mz)
    spsall_pmz <- list(sps_all, pre_mz)
    return(spsall_pmz)
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
    
    ####-------------------------------------------------------------
    #### Dereplication with GNPS ----
    ####-------------------------------------------------------------
    
    if (db == "GNPS"){
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(x, spec = "sps_all", ppmx = 15)
        
        #### GNPS spec with pre_mz
        gnps_with_mz <- spec2_Processing(x, spec = "gnps", ppmx = 15)
        
        
        if (length(sps) > 1 && length(gnps_with_mz) >1){
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            gnps_best_match <- gnps_with_mz[idx[2]]
            df_peaklists <- peakdf(gnps_best_match, sps[idx[1]])
            
            if (!(is.null(df_peaklists))){
                
                print("more spectra and more gnps spectra")
                
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
                print("NO - more spectra and more gnps spectra")
                
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
                print("one spectra and more gnps spectra")
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
                print("NO - one spectra and more gnps spectra")
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
                print("more spectra and one gnps spectra")
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
                print("NO - more spectra and one gnps spectra")
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
                print("one spectra and one gnps spectra")
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
                print("NO - one spectra and one gnps spectra")
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
            print("NO Results")
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
    
    ####-------------------------------------------------------------
    #### Dereplication with HMDB ----
    ####-------------------------------------------------------------
    
    if (db == "HMDB"){
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(x, spec = "sps_all", ppmx = NULL)
        
        #### HMDB spec with pre_mz
        hmdb_with_mz <- spec2_Processing(x, spec = "hmdb", ppmx = 15)
        
        if (length(sps) > 1 && length(hmdb_with_mz) > 1){
            #' Compare experimental spectra against HMDB
            res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)
            #' obtain HMDB spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            hmdb_best_match <- hmdb_with_mz[idx[2]]
            df_peaklists <- peakdf(hmdb_best_match, sps[idx[1]])
            
            if (!(is.null(df_peaklists))){
                print("more sps and more hmdb_with_mz")
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
                HMDBcompoundID <- hmdb_best_match$compound_id
                HMDBmirrorSpec <- str_replace(name_plotmirror, input_dir, "./")
                Source <- "HMDB"
            }
            else{
                print("NO - more sps and more hmdb_with_mz")
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
            
            #' obtain HMDB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            hmdb_best_match <- hmdb_with_mz[gx]

            df_peaklists <- peakdf(hmdb_best_match, sps)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                print("one sps and more hmdb_with_mz")
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
                print("NO one sps and more hmdb_with_mz")
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
            #' obtain hmdb spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(hmdb_with_mz, sps[gx])
            hmdb_best_match <- hmdb_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                print("more sps and one hmdb_with_mz")
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
                print("NO more sps and one hmdb_with_mz")
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
                print("one sps and one hmdb_with_mz")
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
                print("NO one sps and one hmdb_with_mz")
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
            print("NO sps and NO hmdb_with_mz")
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
    
    ####-------------------------------------------------------------
    #### Dereplication with MassBank ----
    ####-------------------------------------------------------------
    
    if (db == "MassBank"){
        
        nx <- 0
        nx <- nx+1
        
        #### input spec with pre_mz
        sps <- spec2_Processing(x, spec = "sps_all", ppmx = NULL)
        
        #### GNPS spec with pre_mz
        mbank_with_mz <- spec2_Processing(x, spec = "mbank", ppmx = 15)
        
        
        if (length(sps) > 1 && length(mbank_with_mz) >1){
            #' Compare experimental spectra against GNPS
            res <- compareSpectra(sps, mbank_with_mz, ppm = 15)
            #' obtain GNPS spectra that matches the most with m/z MS2 spectra
            idx <- which(res == max(res), arr.ind = TRUE)
            mbank_best_match <- mbank_with_mz[idx[2]]
            df_peaklists <- peakdf(mbank_best_match, sps[idx[1]])
            
            if (!(is.null(df_peaklists))){
                print("more spectra and more mbank spectra")
                dir_name <- paste(result_dir, "/spectral_dereplication/MassBank/", sep = "")
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
                print("NO - more spectra and more mbank spectra")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                 
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
            #' obtain MassBank spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            mbank_best_match <- mbank_with_mz[gx]

            df_peaklists <- peakdf(mbank_best_match, sps)

            
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                print("one spectra and more mbank spectra")
                dir_name <- paste(result_dir, "/spectral_dereplication/MassBank/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
                print("NO - one spectra and more mbank spectra")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                 
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
            #' obtain MB spectra that matches the most with m/z MS2 spectra
            gx <- which(res == max(res))
            gx <- gx[1]
            sps <- sps[gx]
            df_peaklists <- peakdf(mbank_with_mz, sps[gx])
            mbank_best_match <- mbank_with_mz
            #' if there are more than 2 peak matching
            if (!(is.null(df_peaklists))){
                print("more spectra and one mbank spectra")
                dir_name <- paste(result_dir, "/spectral_dereplication/MassBank/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
                print("NO - more spectra and one mbank spectra")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                 
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
            mbank_best_match <- mbank_with_mz
            df_peaklists <- peakdf(mbank_best_match, sps)
            if (!(is.null(df_peaklists))){
                print("one spectra and one mbank spectra")
                dir_name <- paste(result_dir, "/spectral_dereplication/MassBank/", sep = "")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
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
                print("NO - one spectra and one mbank spectra")
                id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
                premz <- x
                rtmin <- max(sps$rtime)
                rtmax <- min(sps$rtime)
                rtmed <- median(sps$rtime, na.rm = TRUE)
                 
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
            print("NO results in MassBank")
            id_X <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                      "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                      "ID", as.character(nx), sep = '')
            premz <- x
            rtmin <- max(sps$rtime)
            rtmax <- min(sps$rtime)
            rtmed <- median(sps$rtime, na.rm = TRUE)
                 
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
        
        df_mbank <- cbind(id_X, premz, rtmin, rtmax, rtmed, MBmax_similarity, MBmzScore, MBintScore, MQMatchingPeaks, 
                         MBTotalPeaks, mQueryTotalPeaks, MBformula, MBinchiKEY, MBspectrumID, MBcompound_name, MBmirrorSpec, 
                         Source)
        return(df_mbank)
    }
    
}

save.image(file = "R_Functions.RData")



#' Extract MS2 Fragment peaks
# This functon returns a dataframe and stores a csv file 
    # the directory for csv file is input_dir + /insilico/MS2DATA.csv
# input is from spec_Processing and result directory for each mzML input file

ms2_peaks <- function(x, result_dir){
    
    ## Define variables
    premz <- c() # stores mz
    rtmin <- c() # stores rtmin
    rtmax <- c() # stores rtmax
    rtmed <- c() # stores calculated median of rtmin and rtmax
    col_eng <- c() # stores collision energy
    pol <- c() # stores polarity
    ms2Peaks <- c() # stores the peak list file directory
    id_X <- c() # creates a unique ID based on mz, rt and also the index 
        #(since the mz and rt can be similar in some cases)
    int <- c() # stores intensity of the MS1 feature
    nx <- 0 # stores number for the ID 
    indeX <- 0 # stores number to name the peaklist files 
    
    # x[[2]] is a list of precursor m/z
    for (i in x[[2]]){
        #mz
        premz <- c(premz, i)
        
        #filter based on pre mz; x[[1]] is preprocessed spectra
        sps <- filterPrecursorMz(x[[1]], i)
        
        #rtmin
        rn <- min(sps$rtime)
        rtmin <- c(rtmin, rn)
        
        #rtmax
        rx <- max(sps$rtime)
        rtmax <- c(rtmax, rx)
        
        #rtmedian
        rtm <- median(sps$rtime, na.rm = TRUE)
        rtmed <- c(rtmed, rtm)
        
        #collision energy
        ce <- max(sps$collisionEnergy)
        col_eng <- c(col_eng, ce)
            
        #polarity
        pl <- max(sps$polarity)
        if (pl == 1){
            px <- 'pos'
            pol <- c(pol, px)
        }
        else {
            px <- 'neg'
            pol <- c(pol, px)
        }
        
        #int 
        ints <- max(sps$precursorIntensity)
        int <- c(int, ints) 
        
        #ids
        nx <- nx+1
        fileR <- paste("M", as.character(round(i, digits = 0)),"R",as.character(round(rtm, digits = 0)), "ID", as.character(nx), sep = '')
        id_X <- c(id_X, fileR)
        
        #peak lists
        # variable for name
        names <- c()
        
        # create a new directory to store all the peak list txt files
        if (!file.exists(paste(result_dir, "/insilico/peakfiles_ms2", sep =""))){
            dir.create(paste(result_dir, "/insilico/peakfiles_ms2", sep =""), recursive = TRUE)
        }
        
        for (j in 1:length(sps)){
            nam <- paste('pk', j, sep = '') ## name of variable
            assign(nam, cbind(mz = unlist(mz(sps[j])),intensity = unlist(intensity(sps[j])))) ## assign name to peaklist
            names <- c(names, nam) ## save names in another variable
        
            ## at the end of each list, extract the peak list via combinePeaks function
            if (j == length(sps)){
                n <- paste(names, collapse = ', ') #paste names at the end
                func <- eval(parse(text = paste('combinePeaks(list(',n,'))', sep = ''))) #write the function and then run it
                indeX <- indeX+1
                Y <- as.character(indeX)# numbering for naming peak lists
                #create separate folder for peaklists files
                fileN <- paste(result_dir, '/insilico/peakfiles_ms2/Peaks_0', Y, '.txt', sep = '')
                write.table(func, fileN, row.names = FALSE, col.names = FALSE)
                fileN1 <- str_replace(fileN, input_dir, "./")
                ms2Peaks <- c(ms2Peaks, fileN1)
            }
        }
    }
    first_list <- data.frame(cbind(id_X, premz, rtmed, int ,col_eng, pol, ms2Peaks))
    return(first_list)
}
# Usage
#spec_pr2 <- ms2_peaks(spec_pr, './MZML/DS_201124_SC_full_PRM_neg_09/')
#spec_pr2

# Extract isotopic peaks for each pre_mz
# The input is x = first_list (from ms2peaks function) and y = camera results 

ms1_peaks <- function(x, y, result_dir, QCfile = TRUE){
    # store the ms1_peak list path here
    ms1Peaks <- c()
    
    if (QCfile){
        
        # create a new directory to store all the peak list txt files
        if (!file.exists(paste(result_dir, "/insilico/peakfiles_ms1", sep =""))){
            dir.create(paste(result_dir, "/insilico/peakfiles_ms1", sep =""), recursive = TRUE)
        }
        
        # read the CAMERA results
        y = read.csv(y)
        
        # for all indices in the ms2 features table
        for (i in 1:nrow(x)){
            #store the indices of CAMERA that have same mz and rt as the ms2 features table
            store_c <- c()
            # for all indices in CAMERA Results
            for (j in 1:nrow(y)){
                # if mz and rt from ms2 features are within the range of CAMERA mz and rt
                if (x[i, 'premz'] <= y[j, "mzmax"]   && y[j, "mzmin"] <= x[i, 'premz'] && x[i, 'rtmed'] <= y[j, "rtmax"] && y[j, "rtmin"] <= x[i, 'rtmed']){
                    store_c <- c(store_c, j)
                }
            }
            # indices with same pre m/z and same rt
            df_y <- y[store_c, ]
            df_y <- as.data.frame(df_y)
            #print(df_y)
            #if there was only one index
            if (nrow(df_y) == 1){
                # if there was no isotope annotation for that one index
                if (is.na(df_y[1, "isotopes"])){
                
                    mz <- df_y[1, "mz"] # save mz
                    int <- df_y[1, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(result_dir, "/insilico/peakfiles_ms1/ms1_peaks_", x[i, 'premz'], "_no_isotopes.txt", sep = "") # save name of the peaklist
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE) # save peak list
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1) # add the path of the peak list to a list
                }
                # if there was an isotope annotation
                else{
                
                    df_x <- y[which(y[, "file_origin"] ==df_y[1, "file_origin"]), ] # extract camera results from one file origin
                    df_x <- df_x[which(df_x[, 'istops'] == df_y[1, 'istops']), ] # extract only certain isotope annotation group
                    mz <- df_x[, "mz"] # save mz
                    int <- df_x[, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(result_dir, "/insilico/peakfiles_ms1/ms1_peaksISOTOPE_", x[i, 'premz'], "_isotopeNum_", df_x[1, "istops"], ".txt", sep = "")
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE)
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1)
                }
            }
            # if there are more indices for df_y
            else if(nrow(df_y) > 1){
                # if all enteries have no isotope annotation
                if(all(is.na(df_y[, 'istops']))){
                
                    df_z <- df_y[which(df_y[,"into"] == max(df_y[,"into"])), ] # extract the ms1 peak with highest intensity
                    mz <- df_z[1, "mz"] # save mz
                    int <- df_z[1, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(result_dir, "/insilico/peakfiles_ms1/ms1_peaks_", x[i, 'premz'], "_no_isotopes.txt", sep = "") # save name of the peaklist
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE) # save peak list
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1) # add the path of the peak list to a list
                }
                # if not all isotope annotations are NA
                else if (!(all(is.na(df_y[, 'istops'])))){
                
                    df_y <- df_y[!is.na(df_y$'istops'),] # Remove the NA isotope annotations
                    df_z <- df_y[which(df_y[,"into"] == max(df_y[,"into"])), ] # Select the MS1 peak with highest intensity
                    df_z1 <- y[which(y[, "file_origin"] == df_z[1, "file_origin"]), ]  # extract camera results from one file origin
                    df_z1 <- df_z1[which(df_z1[, 'istops'] == df_z[1, 'istops']), ] # extract only certain isotope annotation group
                    mz <- df_z1[, "mz"] # save mz
                    int <- df_z1[, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(result_dir, "/insilico/peakfiles_ms1/ms1_peaksISOTOPE_", x[i, 'premz'], "_isotopeNum_", df_z1[1, 'istops'],".txt", sep = "") # save name of the peaklist
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE) # save peak list
                    name_file1 <- str_replace(name_file, input_dir, "./")
                    ms1Peaks <- c(ms1Peaks, name_file1) # add the path of the peak list to a list
                }
            }
            else if (nrow(df_y)==0){
                ms1Peaks <- c(ms1Peaks, 'no ms1 peaks in QC')
            }
        }
        second_list <- data.frame(cbind(x, ms1Peaks))
        write.csv(second_list, file = paste(result_dir,'/insilico/MS1DATA.csv', sep = ""))
        return(second_list)
    }
    else{
        ms1Peaks <- c(ms1Peaks, 'no ms1 peaks in QC')
        second_list <- data.frame(cbind(x, ms1Peaks))
        write.csv(second_list, file = paste(result_dir,'/insilico/MS1DATA.csv', sep = ""))
        return(second_list)
    }
    
}

# Usage
#ms1p <- ms1_peaks(x = spec_pr2, y = './MZML/QC/NEW/Combined_Camera_neg.csv', result_dir = './MZML/DS_201124_SC_full_PRM_neg_09', QC = TRUE)
#ms1p

sirius_param <- function(x, result_dir){
    
    if (!file.exists(paste(result_dir, "/insilico/SIRIUS", sep = ""))){
        dir.create(paste(result_dir, "/insilico/SIRIUS", sep = ""), recursive = TRUE) ##create folder
    }
    isotopes <- c() #NA or isotope group number
    sirius_param_file <- c() #input for SIRIUS
    outputNames <- c() #output for SIRIUS with all db
    outputNamesSL <- c() #output for SIRIUS with suspect list as db

    parameter_file <- c()
    par <- 0
    
    a <-0 # counting
    y <- 0 # counting
    z <- 0 # counting
    for (i in 1:nrow(x)){
        
        par <- par+1
        para <- as.character(par) # for numbering
        
        #no MS1 PEAKS and no ISOTOPES
        
        if (x[i, "ms1Peaks"] == 'no ms1 peaks in QC'){

            #INPUT FILE NAME
            fileR <- paste(paste(result_dir, "/insilico/SIRIUS/", sep = ""), para, "_NA_iso_NA_MS1p_", x[i, "premz"], "_SIRIUS_param.ms", sep = "")
            
            sirius_param_file <- c(sirius_param_file, fileR)
            #ISOTOPE Information
            isotopes <- c(isotopes, NA)
            #OUTPUT
            fileSR <- paste(str_remove(fileR, ".ms"),'.json', sep = '')
            fileSRS <- paste(str_remove(fileR, ".ms"), 'SList.json', sep = '')
            outputNames <- c(outputNames, fileSR)
            outputNamesSL <- c(outputNamesSL, fileSRS)
            file.create(fileR, recursive = TRUE)
            file.conn <- file(fileR)
            open(file.conn, open = "at")
            
            #compound
            writeLines(paste(">compound", x[i,"id_X"], sep=" "),con=file.conn)
            #parentmass
            writeLines(paste(">parentmass", x[i,"premz"], sep=" "),con=file.conn)
            ##charge
            if (x[i,"pol"] == "pos"){
                writeLines(paste(">charge", "+1" ,sep=" "),con=file.conn)
            }
            else{
                writeLines(paste(">charge", "-1" ,sep=" "),con=file.conn)
            }
            #rt
            writeLines(paste(">rt", paste(x[i,"rtmed"], "s", sep =''), sep=" "),con=file.conn)
            
            #ms1
            writeLines(">ms1",con=file.conn)
            writeLines(paste(x[i,"premz"], x[i,"int"] ,sep=" "),con=file.conn)
            
            #ms2
            writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            
            peak<- read.table(x[i,"ms2Peaks"])
            for (k in 1:length(peak[,1])){
                writeLines(paste(as.character(peak[k,1]),as.character(peak[k,2]), sep =" "), con=file.conn) 
            }
            close(file.conn)
            parameter_file <- c(parameter_file,file.conn)
        }
        
        # MS1 PEAKS and no ISOTOPES
        
        else if (grepl("_no_isotopes.txt", x[i, "ms1Peaks"], fixed=TRUE)){

            #INPUT FILE NAME
            fileR <- paste(paste(result_dir, "/insilico/SIRIUS/", sep = ""), para, "_NA_iso_MS1p_", x[i, "premz"], "_SIRIUS_param.ms", sep = "")
            sirius_param_file <- c(sirius_param_file, fileR)
            #ISOTOPE Information
            isotopes <- c(isotopes, NA)
            #OUTPUT
            fileSR <- paste(str_remove(fileR, ".ms"),'.json', sep = '')
            fileSRS <- paste(str_remove(fileR, ".ms"), 'SList.json', sep = '')
            outputNames <- c(outputNames, fileSR)
            outputNamesSL <- c(outputNamesSL, fileSRS)
            file.create(fileR, recursive = TRUE)
            file.conn <- file(fileR)
            open(file.conn, open = "at")
            
            #compound
            writeLines(paste(">compound", x[i,"id_X"], sep=" "),con=file.conn)
            #parentmass
            writeLines(paste(">parentmass", x[i,"premz"], sep=" "),con=file.conn)
            ##charge
            if (x[i,"pol"] == "pos"){
                writeLines(paste(">charge", "+1" ,sep=" "),con=file.conn)
            }
            else{
                writeLines(paste(">charge", "-1" ,sep=" "),con=file.conn)
            }
            #rt
            writeLines(paste(">rt", paste(x[i,"rtmed"], "s", sep =''), sep=" "),con=file.conn)
            
            #ms1
            writeLines(">ms1",con=file.conn)
            peakms1<- read.table(x[i,"ms1Peaks"])
            for (l in 1:length(peakms1[,1])){
                writeLines(paste(as.character(peakms1[l,1]),as.character(peakms1[l,2]), sep =" "), con=file.conn) 
            }
            
            #ms2
            writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            
            peakms2<- read.table(x[i,"ms2Peaks"])
            for (k in 1:length(peakms2[,1])){
                writeLines(paste(as.character(peakms2[k,1]),as.character(peakms2[k,2]), sep =" "), con=file.conn) 
            }
            
            close(file.conn)
            parameter_file <- c(parameter_file,file.conn)
        }
        
        # MS1 PEAKS and ISOTOPES
        
        else if (grepl("_isotopeNum_", x[i, "ms1Peaks"], fixed=TRUE)){

            #INPUT FILE NAME
            fileR <- paste(paste(result_dir, "/insilico/SIRIUS/", sep = ""), para, "_iso_MS1p_", as.character(x[i, "premz"]), "_SIRIUS_param.ms", sep = "")
            sirius_param_file <- c(sirius_param_file, fileR)
            #ISOTOPE Information
            isotopes <- c(isotopes, "present")
            #OUTPUT
            fileSR <- paste(str_remove(fileR, ".ms"),'.json', sep = '')
            fileSRS <- paste(str_remove(fileR, ".ms"), 'SList.json', sep = '')
            outputNames <- c(outputNames, fileSR)
            outputNamesSL <- c(outputNamesSL, fileSRS)
            file.create(fileR, recursive = TRUE)
            file.conn <- file(fileR)
            open(file.conn, open = "at")
            
             #compound
            writeLines(paste(">compound", x[i,"id_X"], sep=" "),con=file.conn)
            #parentmass
            writeLines(paste(">parentmass", x[i,"premz"], sep=" "),con=file.conn)
            ##charge
            if (x[i,"pol"] == "pos"){
                writeLines(paste(">charge", "+1" ,sep=" "),con=file.conn)
            }
            else{
                writeLines(paste(">charge", "-1" ,sep=" "),con=file.conn)
            }
            #rt
            writeLines(paste(">rt", paste(x[i,"rtmed"], "s", sep =''), sep=" "),con=file.conn)
            
            
            #ms1
            writeLines(">ms1",con=file.conn)
            peakms1<- read.table(x[i,"ms1Peaks"])
            for (l in 1:length(peakms1[,1])){
                writeLines(paste(as.character(peakms1[l,1]),as.character(peakms1[l,2]), sep =" "), con=file.conn) 
            }
            
            #ms2
            writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            
            peakms2<- read.table(x[i,"ms2Peaks"])
            for (k in 1:length(peakms2[,1])){
                writeLines(paste(as.character(peakms2[k,1]),as.character(peakms2[k,2]), sep =" "), con=file.conn) 
            }
            
            close(file.conn)
            parameter_file <- c(parameter_file,file.conn)
            
        }
    }
    return(data.frame(cbind(sirius_param_file, outputNames, outputNamesSL, isotopes)))
}


#sirius_param_files <- sirius_param(ms1p, result_dir = './MZML/DS_201124_SC_full_PRM_neg_09')

save.image(file = "R_Functions.RData")







save.image(file = "R_Functions.RData")






