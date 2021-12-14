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
    # x is mzML files
    for (i in 1:length(mzml_files)){
        # remove .mzML to extract just the names
        name_mzml <- str_remove(as.character(mzml_files[i]), ".mzML")      
        #' for each file a subdirectory is created to store all results in that, add working directory
        if (!file.exists(name_mzml)){
            dir.create(name_mzml) ##create folder
        }
        ResultFileNames<- c(ResultFileNames, name_mzml)
    }
    input_table <- cbind(mzml_files, ResultFileNames)
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

spec_dereplication<- function(p_mz, db){
    index <- c() ## index
    mz_list <- c() ## m/z
    
    rtmin <- c()
    rtmax <- c()
    
    GNPSmax_similarity <- c() ## GNPS score
    GNPSmzScore <- c()
    GNPSintScore <- c()
    
    GQMatchingPeaks <- c() ##no. of matching peaks
    GNPSTotalPeaks <- c() ##no. of query peaks
    gQueryTotalPeaks<- c()
    
    GNPSSMILES <- c()
    GNPSspectrumID <- c()
    GNPScompound_name <- c() ## Compound name
    
    GNPSdir <- c() ## Result directory
    Source <- c()
    
    if (db = "GNPS"){
        
    }
    
}

    #' Define empty variables for Spectral Screening Results table
    index <- c() ## index
    mz_list <- c() ## m/z
    
    rtmin <- c()
    rtmax <- c()
    
    GNPSmax_similarity <- c() ## GNPS score
    GNPSmzScore <- c()
    GNPSintScore <- c()
    
    GQMatchingPeaks <- c() ##no. of matching peaks
    GNPSTotalPeaks <- c() ##no. of query peaks
    gQueryTotalPeaks<- c()
    
    GNPSSMILES <- c()
    GNPSspectrumID <- c()
    GNPScompound_name <- c() ## Compound name
    
    GNPSresult_dir <- c() ## Result directory
    Results <- c()


    x <- 0 ## to index
    
    ###################################################################
    
    ## SECOND LOOP -- 1.1
    
    #' for all m/z(s) in the input mzML file
    for(a in pre_mz){
        
        ## index
        x <- x+1
        index <- c(index, x)

        ## m/z_list
        mzs <- a
        mz_list <- c(mz_list, mzs)

        ##----------------------------------------------------------------
        ## Spectral Screening with GNPS -
        ##----------------------------------------------------------------
        
        sps <- spec2_Processing(a, spec = "sps_all", ppmx = NULL)
        ## rtmin
        rn <- min(sps$rtime)
        rtmin <- c(rtmin, rn)
        
        ## rtmax
        rx <- max(sps$rtime)
        rtmax <- c(rtmax, rx)
        
        #### GNPS
        gnps_with_mz <- spec2_Processing(a, spec = "gnps", ppmx = 15)
        
        #' Compare experimental spectra against GNPS
        res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)
        
        if(max(res) == 0 || is.nan(max(res))){
            lists <- list(GNPSmax_similarity <- c(GNPSmax_similarity, NA), GNPSintScore <- c(GNPSintScore, NA)
                          , GNPSmzScore <- c(GNPSmzScore, NA), GNPSresult_dir <- c(GNPSresult_dir, NA), 
                          GQMatchingPeaks <- c(GQMatchingPeaks, NA), gQueryTotalPeaks<- c(gQueryTotalPeaks, NA)
                          , GNPSTotalPeaks <- c(GNPSTotalPeaks, NA), GNPScompound_name <- c(GNPScompound_name, NA)
                          , GNPSSMILES <- c(GNPSSMILES, NA), GNPSspectrumID <- c(GNPSspectrumID, NA))
            Results <- c(Results, "no GNPS spectral matching")
        }
        else{ 
            if (length(sps) > 1 && length(gnps_with_mz) >1){
                #' obtain GNPS spectra that matches the most with m/z MS2 spectra
                idx <- which(res == max(res), arr.ind = TRUE)
                gnps_best_match <- gnps_with_mz[idx[2]]
                
                df_peaklists <- peakdf(gnps_best_match, sps[idx[1]])
            
                #' if there are more than 2 peak matching
                if (!(is.null(df_peaklists))){
                
                    #' Identify the best-matching pair
                    mx <- max(res)
                    GNPSmax_similarity <- c(GNPSmax_similarity, mx)
                
                    int_score <- mean(1-(df_peaklists[,"diff"]/100))
                    GNPSintScore <- c(GNPSintScore, int_score)
                    
                    mz_score <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps[idx[1]])[[1]]))
                    GNPSmzScore <- c(GNPSmzScore, mz_score)
                
                
                    Mpeaks <- nrow(df_peaklists)
                    GQMatchingPeaks <- c(GQMatchingPeaks, Mpeaks) ##no. of matching peaks
                
                    GnpsPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GnpsPeaks) ##no. of query peaks
                
                    QueryPeaks <- nrow(peaksData(sps[idx[1]])[[1]]) #total no. of query peaks
                    gQueryTotalPeaks<- c(gQueryTotalPeaks, QueryPeaks)
                    
                    Results <- c(Results, "GNPS")
                    
                    #' compound name
                    comp <- gnps_best_match$NAME
                    GNPScompound_name <- c(GNPScompound_name, comp)
                    id_gnps<- gnps_best_match$SPECTRUMID
                    GNPSspectrumID <- c(GNPSspectrumID, id_gnps)
                    smiles_gnps <- gnps_best_match$SMILES
                    GNPSSMILES<- c(GNPSSMILES, smiles_gnps)
                
                    #' GNPS results directory
                    dir_name <- paste(GNPSdirName, "/mz_", mzs, "/", sep ="")
                    dir.create(dir_name)
                    GNPSresult_dir <- c(GNPSresult_dir, dir_name)
                
                    #'  mgf file

                    #name_mgf <- paste(dir_name, mz,"_", mbank_best_match$compound_name, ".mgf", sep ="")
                    #export(sps, backend=MsBackendMgf(), file = name_mgf)

                    #' plotMirror
                    name_plotmirror <- paste(dir_name, mzs,"_spectra_vs_", gnps_best_match$SPECTRUMID, "_spectra.pdf", sep ="")
                    pdf(name_plotmirror)
                    plotSpectraMirror(sps[idx[1]], gnps_with_mz[idx[2]], tolerance = 0.2,
                                      labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                      labelSrt = -30)
                    grid()
                    dev.off() 
        
                } # if the nrow is 1 
                else{
                    lists <- list(GNPSmax_similarity <- c(GNPSmax_similarity, NA), GNPSintScore <- c(GNPSintScore, NA)
                          , GNPSmzScore <- c(GNPSmzScore, NA), GNPSresult_dir <- c(GNPSresult_dir, NA), 
                          GQMatchingPeaks <- c(GQMatchingPeaks, NA), gQueryTotalPeaks<- c(gQueryTotalPeaks, NA)
                          , GNPSTotalPeaks <- c(GNPSTotalPeaks, NA), GNPScompound_name <- c(GNPScompound_name, NA)
                          , GNPSSMILES <- c(GNPSSMILES, NA), GNPSspectrumID <- c(GNPSspectrumID, NA))
                    Results <- c(Results, "no matching peaks")
                }

        } else if (length(sps) == 1 && length(gnps_with_mz) >1){
                #' obtain GNPS spectra that matches the most with m/z MS2 spectra
                gx <- which(res == max(res))
                gx <- gx[1]
                gnps_best_match <- gnps_with_mz[gx]

                df_peaklists <- peakdf(gnps_best_match, sps)

            
                #' if there are more than 2 peak matching
                if (!(is.null(df_peaklists))){
                
                    #' Identify the best-matching pair
                    mx <- max(res)
                    GNPSmax_similarity <- c(GNPSmax_similarity, mx)
                
                    int_score <- mean(1-(df_peaklists[,"diff"]/100))
                    GNPSintScore <- c(GNPSintScore, int_score)
                    
                    mz_score <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    GNPSmzScore <- c(GNPSmzScore, mz_score)
                
                
                    Mpeaks <- nrow(df_peaklists)
                    GQMatchingPeaks <- c(GQMatchingPeaks, Mpeaks) ##no. of matching peaks
                    
                    GnpsPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GnpsPeaks) ##no. of query peaks
                
                    QueryPeaks <- nrow(peaksData(sps)[[1]])
                    gQueryTotalPeaks<- c(gQueryTotalPeaks, QueryPeaks)
                    
                    Results <- c(Results, "GNPS")
                
                    #' compound name
                    comp <- gnps_best_match$NAME
                    GNPScompound_name <- c(GNPScompound_name, comp)
                    id_gnps<- gnps_best_match$SPECTRUMID
                    GNPSspectrumID <- c(GNPSspectrumID, id_gnps)
                    smiles_gnps <- gnps_best_match$SMILES
                    GNPSSMILES<- c(GNPSSMILES, smiles_gnps)
                
                    #' GNPS results directory
                    dir_name <- paste(GNPSdirName, "/mz_", mzs, "/", sep ="")
                    dir.create(dir_name)
                    GNPSresult_dir <- c(GNPSresult_dir, dir_name)

                    #' plotMirror
                    name_plotmirror <- paste(dir_name, mzs,"_spectra_vs_", gnps_best_match$SPECTRUMID, "_spectra.pdf", sep ="")
                    pdf(name_plotmirror)
                    plotSpectraMirror(sps, gnps_best_match, tolerance = 0.2,
                                      labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                      labelSrt = -30)
                    grid()
                    dev.off() 
                } # if the nrow is 1 
                else{
                    lists <- list(GNPSmax_similarity <- c(GNPSmax_similarity, NA), GNPSintScore <- c(GNPSintScore, NA)
                          , GNPSmzScore <- c(GNPSmzScore, NA), GNPSresult_dir <- c(GNPSresult_dir, NA), 
                          GQMatchingPeaks <- c(GQMatchingPeaks, NA), gQueryTotalPeaks<- c(gQueryTotalPeaks, NA)
                          , GNPSTotalPeaks <- c(GNPSTotalPeaks, NA), GNPScompound_name <- c(GNPScompound_name, NA)
                          , GNPSSMILES <- c(GNPSSMILES, NA), GNPSspectrumID <- c(GNPSspectrumID, NA))
                    Results <- c(Results, "no matching peaks")
                }
            } else if (length(sps) > 1 && length(gnps_with_mz) == 1){
                #' obtain MB spectra that matches the most with m/z MS2 spectra
                gx <- which(res == max(res))
                gx <- gx[1]
                sps <- sps[gx]
                df_peaklists <- peakdf(gnps_with_mz, sps[gx])
                gnps_best_match <- gnps_with_mz

                #' if there are more than 2 peak matching
                if (!(is.null(df_peaklists))){
                
                    #' Identify the best-matching pair
                    mx <- max(res)
                    GNPSmax_similarity <- c(GNPSmax_similarity, mx)
                
                    int_score <- mean(1-(df_peaklists[,"diff"]/100))
                    GNPSintScore <- c(GNPSintScore, int_score)
                    
                    mz_score <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_best_match)[[1]])+nrow(peaksData(sps)[[1]]))
                    GNPSmzScore <- c(GNPSmzScore, mz_score)
                
                    Mpeaks <- nrow(df_peaklists)
                    GQMatchingPeaks <- c(GQMatchingPeaks, Mpeaks) ##no. of matching peaks
                    
                    GnpsPeaks <- nrow(peaksData(gnps_best_match)[[1]])
                    GNPSTotalPeaks <- c(GNPSTotalPeaks, GnpsPeaks) ##no. of query peaks
                
                    QueryPeaks <- nrow(peaksData(sps)[[1]])#total no. of query peaks
                    gQueryTotalPeaks<- c(gQueryTotalPeaks, QueryPeaks)
                    
                    Results <- c(Results, "GNPS")
                
                    #' compound name
                    comp <- gnps_best_match$NAME
                    GNPScompound_name <- c(GNPScompound_name, comp)
                    id_gnps<- gnps_best_match$SPECTRUMID
                    GNPSspectrumID <- c(GNPSspectrumID, id_gnps)
                    smiles_gnps <- gnps_best_match$SMILES
                    GNPSSMILES<- c(GNPSSMILES, smiles_gnps)
                
                    #' GNPS results directory
                    dir_name <- paste(GNPSdirName, "/mz_", mzs, "/", sep ="")
                    dir.create(dir_name)
                    GNPSresult_dir <- c(GNPSresult_dir, dir_name)

                    #' plotMirror
                    name_plotmirror <- paste(dir_name, mzs,"_spectra_vs_", gnps_best_match$SPECTRUMID, "_spectra.pdf", sep ="")
                    pdf(name_plotmirror)
                    plotSpectraMirror(sps, gnps_best_match, tolerance = 0.2,
                                      labels = label_fun, labelPos = 2, labelOffset = 0.2,
                                      labelSrt = -30)
                    grid()
                    dev.off() 
                }  
                else{
                    lists <- list(GNPSmax_similarity <- c(GNPSmax_similarity, NA), GNPSintScore <- c(GNPSintScore, NA)
                          , GNPSmzScore <- c(GNPSmzScore, NA), GNPSresult_dir <- c(GNPSresult_dir, NA), 
                          GQMatchingPeaks <- c(GQMatchingPeaks, NA), gQueryTotalPeaks<- c(gQueryTotalPeaks, NA)
                          , GNPSTotalPeaks <- c(GNPSTotalPeaks, NA), GNPScompound_name <- c(GNPScompound_name, NA)
                          , GNPSSMILES <- c(GNPSSMILES, NA), GNPSspectrumID <- c(GNPSspectrumID, NA))
                    Results <- c(Results, "no matching peaks")
                }
            } else {
                lists <- list(GNPSmax_similarity <- c(GNPSmax_similarity, NA), GNPSintScore <- c(GNPSintScore, NA)
                          , GNPSmzScore <- c(GNPSmzScore, NA), GNPSresult_dir <- c(GNPSresult_dir, NA), 
                          GQMatchingPeaks <- c(GQMatchingPeaks, NA), gQueryTotalPeaks<- c(gQueryTotalPeaks, NA)
                          , GNPSTotalPeaks <- c(GNPSTotalPeaks, NA), GNPScompound_name <- c(GNPScompound_name, NA)
                          , GNPSSMILES <- c(GNPSSMILES, NA), GNPSspectrumID <- c(GNPSspectrumID, NA))
                    Results <- c(Results, "no spectral matching")
            }
        }
    }
    gnps_df <- cbind(index, mz_list, rtmin, rtmax, GNPSmax_similarity, GNPSmzScore, GNPSintScore, GQMatchingPeaks, gQueryTotalPeaks, GNPSTotalPeaks, GNPScompound_name, GNPSSMILES, GNPSspectrumID, GNPSresult_dir, Results)
    gnps_df <- data.frame(gnps_df)
    gnps_results_csv <- paste(sub_directory_mzmlfiles[i], "/gnps_csv.csv", sep ="")
    write.csv(gnps_df, gnps_results_csv)




















save.image(file = "R_Functions.RData")




