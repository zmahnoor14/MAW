download_specDB <- function(input_dir, db = "all", error = TRUE){
    # Track Time 
    start_time <- Sys.time()

    # only input available as of now
    databases <- 'gnps, hmdb, mbank, all'
    
    # creat a summary file, open and store timings of download and version if possible
    summaryFile <- paste(input_dir, "summaryFile.txt", sep = "")
    file.create(summaryFile, recursive = TRUE)
    file.conn <- file(summaryFile)
    open(file.conn, open = "at")
            
    # gnps
    if (db == "all" || db =="gnps"){
        
        print("GNPS WORKS")
        
        # Download file
        system(paste("wget -P", 
                     input_dir,
                     "https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.mgf",
                     sep =  " "))
        
        # load the spectra into MsBackendMgf
        gnps <- Spectra(paste(input_dir, "ALL_GNPS.mgf", sep = ''), source = MsBackendMgf())
        save(gnps, file = paste(input_dir,"gnps.rda", sep = ""))
        
        # delete the database in its format to free up space
        system(paste("rm", (paste(input_dir, "ALL_GNPS.mgf", sep = '')), sep = " "))
        
        writeLines(paste("GNPS saved at", Sys.time(), sep=" "),con=file.conn)
        
    }
    # hmdb
    if (db == "all" || db =="hmdb"){
        
        print("HMDB WORKS")
        
        
        
        ####### Version Control ######
        
        # extract HMDB Current version
        html <- read_html("https://hmdb.ca/downloads")
        strings <- html%>% html_elements("a") %>% html_text2()
        ls <- unique(strings)
        hmdb_curr_ver <- c()
        for (i in ls){
            if (grepl("Current", i)){
            hmdb_curr_ver<- c(i, hmdb_curr_ver)
            }
        }
        
        
        
        
        ####### Download and unzip ######
        
        #Download file predicted MSMS spectra
        system(paste("wget",
                     "https://hmdb.ca/system/downloads/current/spectral_data/spectra_xml/hmdb_predicted_msms_spectra.zip",
                     sep = " "))
        # unzip
        system(paste("unzip", "hmdb_predicted_msms_spectra.zip", "-d",  paste(input_dir, "hmdb_predicted_msms_spectra", sep = ""), sep = " "))
    
        
        #Download file experimental MSMS spectra
        system(paste("wget",
                     "https://hmdb.ca/system/downloads/current/spectral_data/spectra_xml/hmdb_experimental_msms_spectra.zip",
                     sep = " "))
        # unzip
        system(paste("unzip", "hmdb_experimental_msms_spectra.zip", "-d", paste(input_dir, "hmdb_experimental_msms_spectra", sep = ""), sep = " "))
        
        
        
        
        ####### Load spectra in MsBackend #######
        
        hmdb_predfiles <- list.files(path = paste(input_dir, "hmdb_predicted_msms_spectra", sep = ''), full.names = TRUE)
        
        hmdb_predicted <- c()
        for (i in hmdb_predfiles){
            # load the spectra into MsBackendHMDB
            hmdb_pred <- Spectra(i, source = MsBackendHmdbXml())
            hmdb_predicted <- c(hmdb_predicted, hmdb_pred)
        }
        
        hmdb_expfiles <- list.files(path = paste(input_dir, "hmdb_experimental_msms_spectra", sep = ''), full.names = TRUE)
        
        hmdb_experimental <- c()
        
        for (j in hmdb_expfiles){
            # load the spectra into MsBackendHMDB
            hmdb_exp <- Spectra(j, source = MsBackendHmdb())
            hmdb_experimental <- c(hmdb_experimental, hmdb_exp)
        }
        
        
        hmdb <- hmdb_predicted + hmdb_experimental
        save(hmdb, file = paste(input_dir,"hmdb.rda", sep = ""))
        
        
        
        
        
        ####### Remove the XML files #######
        
        # delete the database in its format to free up space
        system(paste("rm -r", (paste(input_dir, "hmdb_predicted_msms_spectra", sep = '')), sep = " "))
        system(paste("rm -r", (paste(input_dir, "hmdb_experimental_msms_spectra", sep = '')), sep = " "))
        
        
        writeLines(paste("HMDB saved at", Sys.time(), "with release version", hmdb_curr_ver, sep=" "),con=file.conn)
    }
    
    #mbank
    if (db == "all" || db =="mbank"){
        
        print("MassBank WORKS")
        
        page <- read_html("https://github.com/MassBank/MassBank-data/releases")
        page %>%
            html_nodes("a") %>%       # find all links
            html_attr("href") %>%     # get the url
            str_subset("MassBank_NIST.msp") -> tmp # find those that have the name MassBank_NIST.msp
        
        #download file
        system(paste("wget ",
                     "https://github.com", tmp[1], 
                     sep =  ""))
        
        mbank <- Spectra(paste(input_dir, "MassBank_NIST.msp", sep = ''), source = MsBackendMsp())
        save(mbank, file = paste(input_dir,"mbankNIST.rda", sep = ""))
        
        # delete the database in its format to free up space
        system(paste("rm", (paste(input_dir, "MassBank_NIST.msp", sep = '')), sep = " "))
        
        # obtain the month and year for the database release to add to summary
        res <- str_match(tmp[1], "download/\\s*(.*?)\\s*/MassBank_NIST")
        
        writeLines(paste("MassBank saved at", Sys.time(), "with release version", res[,2], sep=" "),con=file.conn)
    }
    
    #wrong input error message
    else if (!grepl(db, databases, fixed = TRUE)){
        stop("Wrong db input. Following inputs apply: gnps, hmdb, mbank or all")
    }
    close(file.conn)
    #download_specDB(input_dir, db)
    end_time <- Sys.time()
    print(end_time - start_time)
}




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

# Usage:
# addProcessing(sps, norm_int)

#' addProcessing is a predefined function in Spectra package

##-----------------------------------------------------------------
## Result Directories 
##-----------------------------------------------------------------

## Specifying a function for creating result directories for each input mzml
# input for the function:
# input directory
ms2_rfilename<- function(input_dir){
    #list_ms2_files <- intersect(list.files(input_dir, pattern = "_PRM_"), list.files(input_dir, pattern = ".mzML"))
    list_ms2_files <- list.files(input_dir, pattern = ".mzML")
    mzml_file <- paste(input_dir, list_ms2_files, sep = "")
    
    #store the result file names to return to this function as output
    mzml_files <- c()
    ResultFileNames <- c()
    File_id <- c()
    nx <- 0
    # x is mzML files
    for (i in 1:length(mzml_file)){
        nx <- nx+1
        # remove .mzML to extract just the names
        mzml_filex <- str_replace(mzml_file[i], input_dir, "./")
        name_mzmls <- str_remove(as.character(mzml_filex), ".mzML")
        #name_mzml <- str_replace(name_mzmls, input_dir, "./")
        #' for each file a subdirectory is created to store all results in that, add working directory
        if (!file.exists(name_mzmls)){
            dir.create(name_mzmls) ##create folder
        }
        ResultFileNames<- c(ResultFileNames, name_mzmls)
        mzml_files <- c(mzml_files, mzml_filex)
        File_id <- c(File_id, paste("file_", nx, sep = ""))
    }
    input_table <- cbind(mzml_files, ResultFileNames, File_id)
    
    write.csv(input_table, paste(input_dir, "input_table.csv", sep = ""))
    return(data.frame(input_table))
}

# usage:
## input directory ##
#input_dir <- paste(getwd(), "/", sep = '')
#ms2_rfilename(input_dir)



##-----------------------------------------------------------------
## Read mzML files and extract precursor m/z(s)
##-----------------------------------------------------------------

#' All spectra in mzML files preprocessing, return two outputs, pre-processed MS2 spectra and all precursor masses
# x is one mzML file
spec_Processing <- function(x, result_dir){
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
    export(sps_all, backend = MsBackendMzR(), file = paste(result_dir, "/processedSpectra.mzML", sep = ""))
    write.table(pre_mz, file = paste(result_dir, "/premz_list.txt", sep = ""), sep = "/t",row.names = FALSE, col.names = FALSE)
    spsall_pmz <- list(sps_all, pre_mz)
    return(spsall_pmz)
}

##-----------------------------------------------------------------
## Pre-process MS2 spectra
##-----------------------------------------------------------------

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
# Usage:
# peakdf(sps[[1]], gnps_best_match, ppmx = 15)

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

# x is one pre_mz, db is GNPS, HMDB, MassBank
spec_dereplication<- function(pre_tbl, proc_mzml, db, result_dir, file_id, input_dir, ppmx, error = TRUE){
    
    ####-------------------------------------------------------------
    #### Dereplication with all or GNPS ----
    ####-------------------------------------------------------------
    
    databases <- 'gnps, hmdb, mbank, all'
    
    sps_all <- Spectra(proc_mzml, backend = MsBackendMzR())
        
    tbl <- read.table(pre_tbl)
    pre_mz <- tbl[[1]]
    
    if (db == "all" || db =="gnps"){
        
        
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
            
            spsrt <- filterPrecursorMz(sps_all, x)
        
            
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
        
        load(file = paste(input_dir,"hmdb.rda", sep = ""))
        
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
            
            spsrt <- filterPrecursorMz(sps_all, x)
        
            
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

        load(file = paste(input_dir,"mbank.rda", sep = ""))
        
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
            
            
            spsrt <- filterPrecursorMz(sps_all, x)
        
            
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
    
    #wrong input error message
    else if (!grepl(db, databases, fixed = TRUE)){
        stop("Wrong db input. Following inputs apply: gnps, hmdb, mbank or all")
    }

}
# Usage
# spec_dereplication<- function(x, db, result_dir, file_id, input_dir, ppmx)


#save.image(file = "R_Functions.RData")

#' Extract MS2 Fragment peaks
# This functon returns a dataframe and stores a csv file 
    # the directory for csv file is input_dir + /insilico/MS2DATA.csv
# input is from spec_Processing and result directory for each mzML input file

ms2_peaks <- function(pre_tbl, proc_mzml, input_dir, result_dir, file_id){
    
    
    sps_all <- Spectra(proc_mzml, backend = MsBackendMzR())
        
    tbl <- read.table(pre_tbl)
    pre_mz <- tbl[[1]]
    
    
    ## Define variables
    premz <- c() # stores mz
    rtmin <- c() # stores rtmin
    rtmax <- c() # stores rtmax
    rtmed <- c() # stores calculated median of rtmin and rtmax
    rtmean <- c() # stores calculated mean of rtmin and rtmax
    col_eng <- c() # stores collision energy
    pol <- c() # stores polarity
    ms2Peaks <- c() # stores the peak list file directory
    id_X <- c() # creates a unique ID based on mz, rt and also the index 
        #(since the mz and rt can be similar in some cases)
    int <- c() # stores intensity of the MS1 feature
    nx <- 0 # stores number for the ID 
    indeX <- 0 # stores number to name the peaklist files 
    
    # pre_mz is a list of precursor m/z
    for (i in pre_mz){
        #mz
        premz <- c(premz, i)
        
        #filter based on pre mz; sps_all is preprocessed spectra
        sps <- filterPrecursorMz(sps_all, i)
        
        #rtmin
        rn <- min(sps$rtime)
        rtmin <- c(rtmin, rn)
        
        #rtmax
        rx <- max(sps$rtime)
        rtmax <- c(rtmax, rx)
        
        #rtmedian
        rtm <- median(sps$rtime, na.rm = TRUE)
        rtmed <- c(rtmed, rtm)
        
        #rtmean
        rtme <- mean(sps$rtime, na.rm = TRUE)
        rtmean <- c(rtmean, rtme)
        
        
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
        id_Xx <- paste(file_id,  "M",  as.character(round(x, digits = 0)), 
                          "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)), 
                          "ID", as.character(nx), sep = '')
        id_X <- c(id_X, id_Xx)
        
        #peak lists
        # variable for name
        names <- c()
        
        # create a new directory to store all the peak list txt files
        dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/peakfiles_ms2", sep =""), "./"), sep = "")
        if (!file.exists(dir_name)){
            dir.create(dir_name, recursive = TRUE)
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
                fileN <- paste(dir_name, '/Peaks_0', Y, '.txt', sep = '')
                write.table(func, fileN, row.names = FALSE, col.names = FALSE)
                fileN1 <- str_replace(fileN, input_dir, "./")
                ms2Peaks <- c(ms2Peaks, fileN1)
            }
        }
    }
    first_list <- data.frame(cbind(id_X, premz, rtmed, rtmean, int ,col_eng, pol, ms2Peaks))
    write.csv(first_list, file = paste(input_dir, str_remove(paste(result_dir,'/insilico/MS2DATA.csv', sep = ""), "./"), sep =""))
    return(first_list)
}
# Usage
#spec_pr2 <- ms2_peaks(spec_pr, './MZML/DS_201124_SC_full_PRM_neg_09/')
#spec_pr2

# for QC files, which have both polarities in one file or the polarities are unknown
# inputs:
#path e.g: “/users/name/project/QC”
#pattern (define this so that the function only catches the file 
#with a certain pattern e.g: “common”, if no pattern, by default 
#function takes “.mzML” as the pattern which is all the .mzML files, 
#so make sure only the QC files are in this directory and not the MS2.mzML files)

cam_funcMode <- function(path, pattern = ".mzML"){
    library("CAMERA")
    # List all files present in QC folder 
    files_QC_N <- list.files(path, pattern = pattern ,full.names=TRUE)
    
    for (i in 1:length(files_QC_N)){
    
        # read each file using Spectra
        sps_all <- Spectra(files_QC_N[i], backend = MsBackendMzR())
        
        
        if (length(unique(sps_all$polarity)) == 1){
            if(unique(sps_all$polarity) == 1){
                #Read the same file with MS1 information; note CAMERA reads xcmsSet object
                xs <- xcmsSet(file = as.character(files_QC_N[i]),
                              profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                              mslevel= 1, progressCallback=NULL, polarity="positive",
                              scanrange = NULL, BPPARAM = bpparam(),
                              stopOnError = TRUE)
                # Create an xsAnnotate object 
                an <- xsAnnotate(xs) 
                # Group based on RT 
                anF <- groupFWHM(an, perfwhm = 0.6)
                # Annotate isotopes 
                anI <- findIsotopes(anF, mzabs = 0.01) 
                # Verify grouping 
                anIC <- groupCorr(anI, cor_eic_th = 0.75)
                #Annotate adducts 
                anFA <- findAdducts(anIC, polarity="positive") 
                #get a feature list
                peaklist <- getPeaklist(anFA)
                # add file_origin information
                peaklist$file_origin <- as.character(files_QC_N[i])
                file_name <- paste(path,"/posCAMERA_Results_", i, ".csv", sep = "")
                # write individual QC files which are in pos mode (later code will combine them)
                write.csv(peaklist, file = file_name)
            }else if (unique(sps_all$polarity) == 0){
                #Read the same file with MS1 information; note CAMERA reads xcmsSet object
                xs <- xcmsSet(file = as.character(files_QC_N[i]),
                              profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                              mslevel= 1, progressCallback=NULL, polarity="negative",
                              scanrange = NULL, BPPARAM = bpparam(),
                              stopOnError = TRUE)
                # Create an xsAnnotate object 
                an <- xsAnnotate(xs) 
                # Group based on RT 
                anF <- groupFWHM(an, perfwhm = 0.6)
                # Annotate isotopes 
                anI <- findIsotopes(anF, mzabs = 0.01) 
                # Verify grouping 
                anIC <- groupCorr(anI, cor_eic_th = 0.75)
                #Annotate adducts 
                anFA <- findAdducts(anIC, polarity="negative") 
                #get a feature list
                peaklist <- getPeaklist(anFA)
                # add file_origin information
                peaklist$file_origin <- as.character(files_QC_N[i])
                file_name <- paste(path, "/negCAMERA_Results_", i, ".csv", sep = "")
                # write individual QC files which are in pos mode (later code will combine them)
                write.csv(peaklist, file = file_name)
            }
        }
    else{
        pos <- sps_all[sps_all$polarity == 1]
        neg <- sps_all[sps_all$polarity == 0]
        file_p <- paste(path, "/QC_280k_pos", i, ".mzML", sep = "")
        file_n <- paste(path, "/QC_280k_neg", i, ".mzML", sep = "")
        
        # create new mzML QC files for pos and neg modes from each common QC file
        export(pos, backend = MsBackendMzR(), file = file_p)
        export(neg, backend = MsBackendMzR(), file = file_n)
        #Read the same file with MS1 information; note CAMERA reads xcmsSet object
        xs <- xcmsSet(file = as.character(file_p),
                        profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                        mslevel= 1, progressCallback=NULL, polarity="positive",
                        scanrange = NULL, BPPARAM = bpparam(),
                        stopOnError = TRUE)
        # Create an xsAnnotate object 
        an <- xsAnnotate(xs) 
        # Group based on RT 
        anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes 
        anI <- findIsotopes(anF, mzabs = 0.01) 
        # Verify grouping 
        anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts 
        anFA <- findAdducts(anIC, polarity="positive") 
        #get a feature list
        peaklist <- getPeaklist(anFA)
        # add file_origin information
        peaklist$file_origin <- as.character(files_QC_N[i])
        file_name <- paste(path, "/posCAMERA_Results_", i, ".csv", sep = "")
        # write individual QC files which are in pos mode (later code will combine them)
        write.csv(peaklist, file = file_name)
        
        #Read the same file with MS1 information; note CAMERA reads xcmsSet object
        xs <- xcmsSet(file = as.character(file_n),
                        profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                        mslevel= 1, progressCallback=NULL, polarity="negative",
                        scanrange = NULL, BPPARAM = bpparam(),
                        stopOnError = TRUE)
        # Create an xsAnnotate object 
        an <- xsAnnotate(xs) 
        # Group based on RT 
        anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes 
        anI <- findIsotopes(anF, mzabs = 0.01) 
        # Verify grouping 
        anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts 
        anFA <- findAdducts(anIC, polarity="negative") 
        #get a feature list
        peaklist <- getPeaklist(anFA)
        # add file_origin information
        peaklist$file_origin <- as.character(files_QC_N[i])
        file_name <- paste(path, "/negCAMERA_Results_", i, ".csv", sep = "")
        # write individual QC files which are in pos mode (later code will combine them)
        write.csv(peaklist, file = file_name)
    }
    
    }
    
    detach("package:CAMERA", unload=TRUE)
}

#Usage: cam_funcMode(path, pattern)

merge_qc<- function(path){
    # combine all QC which are in positive mode
    df_pos <- list.files(path, pattern = "posCAMERA_Results_", full.names = TRUE) %>% 
        lapply(read_csv) %>% 
        bind_rows
    # remove any duplicated rows
    df_pos <- as.data.frame(df_pos[!duplicated(df_pos), ])

    #extract isotope column numbers, the numbers represent the group of isotope
    nm_p <- regmatches(df_pos[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",df_pos[, "isotopes"]))

    # for all the numbers, extract only first number, since it is the group number, 
    # second number can be charge
    for (i in 1:length(nm_p)){
        y <- as.numeric(unlist(nm_p[i]))
        df_pos[i,'istops'] = y[1]
    }

    # write csv for the combined_camera_pos results
    write.csv(df_pos, paste(path, "/Combined_Camera_pos.csv", sep = ""))
    
    # combine all QC which are in negative mode
    df_neg <- list.files(path, pattern = "negCAMERA_Results_", full.names = TRUE) %>% 
        lapply(read_csv) %>% 
        bind_rows
    # remove any duplicated rows based on mz
    df_neg <- as.data.frame(df_neg[!duplicated(df_neg), ])

    #extract isotope column numbers, the numbers represent the group of isotope
    nm_n <- regmatches(df_neg[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",df_neg[, "isotopes"]))

    # for all the numbers, extract only first number, since it is the group number, 
    # second number can be charge
    for (i in 1:length(nm_n)){
        y <- as.numeric(unlist(nm_n[i]))
        df_neg[i,'istops'] = y[1]
    }
    # write csv for the combined_camera_neg results
    write.csv(df_neg, paste(path, "/Combined_Camera_neg.csv", sep = ""))
}
# Usage: merge_qc(path)

cam_func <- function(path, f, mode = "pos"){
    library("CAMERA")
    fl <- paste(path, f, sep ="")
    if(mode == "pos"){
        xs <- xcmsSet(file = fl,profmethod = "bin", 
              profparam = list(), lockMassFreq=FALSE,
              mslevel= 1, progressCallback=NULL, polarity="positive",
              scanrange = NULL, BPPARAM = bpparam(),stopOnError = TRUE)
        # Create an xsAnnotate object 
        an <- xsAnnotate(xs) 
        # Group based on RT 
        anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes 
        anI <- findIsotopes(anF, mzabs = 0.01) 
        # Verify grouping 
        anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts 
        anFA <- findAdducts(anIC, polarity="positive") 
        peaklist <- getPeaklist(anFA)
        peaklist$file_origin <- fl

        #extract isotope column numbers, the numbers represent the group of isotope
        nm_po <- regmatches(peaklist[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",peaklist[, "isotopes"]))
        # for all the numbers in v, extract only first number, since it is the group number, 
        # second number can be charge
        for (i in 1:length(nm_po)){
            y <- as.numeric(unlist(nm_po[i]))
            peaklist[i,'istops'] = y[1]
        }
        name <- str_remove(f, ".mzML")
        write.csv(peaklist, file = paste(input_dir, "/QC/posCAMERAResults_", name,".csv", sep = ""))
    }
    else if(mode == "neg"){
        xs <- xcmsSet(file = fl,profmethod = "bin", 
              profparam = list(), lockMassFreq=FALSE,
              mslevel= 1, progressCallback=NULL, polarity="negative",
              scanrange = NULL, BPPARAM = bpparam(),stopOnError = TRUE)
        # Create an xsAnnotate object 
        an <- xsAnnotate(xs) 
        # Group based on RT 
        anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes 
        anI <- findIsotopes(anF, mzabs = 0.01) 
        # Verify grouping 
        anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts 
        anFA <- findAdducts(anIC, polarity="negative") 
        peaklist <- getPeaklist(anFA)
        peaklist$file_origin <- fl

        #extract isotope column numbers, the numbers represent the group of isotope
        nm_ne <- regmatches(peaklist[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",peaklist[, "isotopes"]))
        # for all the numbers in v, extract only first number, since it is the group number, 
        # second number can be charge
        for (i in 1:length(nm_ne)){
            y <- as.numeric(unlist(nm_ne[i]))
            peaklist[i,'istops'] = y[1]
        }
        name <- str_remove(f, ".mzML")
        write.csv(peaklist, file = paste(input_dir, "QC/negCAMERAResults_", name,".csv", sep = ""))
    }
    detach("package:CAMERA", unload=TRUE)
    
}
# Usage: cam_func(f, mode = "pos")



#addQC_input_table <- function(path, pattern){
#    
#}



# Extract isotopic peaks for each pre_mz
# The input is x = first_list (from ms2peaks function) and y = camera results 

ms1_peaks <- function(x, y, result_dir, input_dir, QCfile = TRUE){
    # store the ms1_peak list path here
    ms1Peaks <- c()
    
    if (QCfile){
        
        dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/peakfiles_ms1", sep =""), "./"), sep = "")
        # create a new directory to store all the peak list txt files
        if (!file.exists(dir_name)){
            dir.create(dir_name, recursive = TRUE)
        }
        
        # read the CAMERA results
        y = read.csv(y)
        x = read.csv(x)
        
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

            #if there was only one index
            if (nrow(df_y) == 1){
                
                # if there was no isotope annotation for that one index
                if (is.na(df_y[1, "istops"])){
                
                    mz <- df_y[1, "mz"] # save mz
                    int <- df_y[1, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(dir_name, "/ms1_peaks_", x[i, 'premz'], "_no_isotopes.txt", sep = "") # save name of the peaklist
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
                    name_file <- paste(dir_name, "/ms1_peaksISOTOPE_", x[i, 'premz'], "_isotopeNum_", df_x[1, "istops"], ".txt", sep = "")
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
                    name_file <- paste(dir_name, "/ms1_peaks_", x[i, 'premz'], "_no_isotopes.txt", sep = "") # save name of the peaklist
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
                    name_file <- paste(dir_name, "/ms1_peaksISOTOPE_", x[i, 'premz'], "_isotopeNum_", df_z1[1, 'istops'],".txt", sep = "") # save name of the peaklist
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
        write.csv(second_list, file = paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA.csv', sep = ""), "./"), sep =""))
        return(second_list)
    }
    
}

# Usage
#ms1p <- ms1_peaks(x = spec_pr2, y = './MZML/QC/NEW/Combined_Camera_neg.csv', result_dir = './MZML/DS_201124_SC_full_PRM_neg_09', QC = TRUE)
#ms1p





# input x is result from either ms1_peaks
# SL is if a suspect list is present

sirius_param <- function(x, result_dir, input_dir, SL = TRUE){
    
    dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/SIRIUS", sep =""), "./"), sep = "")
    if (!file.exists(dir_name)){
        dir.create(dir_name, recursive = TRUE) ##create folder
    }
    isotopes <- c() #NA or isotope group number
    sirius_param_file <- c() #input for SIRIUS
    outputNames <- c() #output for SIRIUS with all db
    outputNamesSL <- c() #output for SIRIUS with suspect list as db

    parameter_file <- c()
    par <- 0
    
    x <- read.csv(x)
    
    a <-0 # counting
    y <- 0 # counting
    z <- 0 # counting
    for (i in 1:nrow(x)){
        
        par <- par+1
        para <- as.character(par) # for numbering
        
        #no MS1 PEAKS and no ISOTOPES
        
        if (x[i, "ms1Peaks"] == 'no ms1 peaks in QC'){

            #INPUT FILE NAME
            fileR <- paste(dir_name, "/" ,para, "_NA_iso_NA_MS1p_", x[i, "premz"], "_SIRIUS_param.ms", sep = "")
            
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
            fileR <- paste(dir_name, "/", para, "_NA_iso_MS1p_", x[i, "premz"], "_SIRIUS_param.ms", sep = "")
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
            fileR <- paste(dir_name, "/", para, "_iso_MS1p_", as.character(x[i, "premz"]), "_SIRIUS_param.ms", sep = "")
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
    if (SL){
        
        write.csv(data.frame(cbind(sirius_param_file, outputNames, outputNamesSL, isotopes)), paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA_SiriusPandSL.csv', sep = ""), "./"), sep =""))
        return(data.frame(cbind(sirius_param_file, outputNames, outputNamesSL, isotopes)))
        
    }
    else{
        
        write.csv(data.frame(cbind(sirius_param_file, outputNames, isotopes)), paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA_SiriusP.csv', sep = ""), "./"), sep =""))
        return(data.frame(cbind(sirius_param_file, outputNames, isotopes)))
        
    }
    
}

# Usage 
# sirius_param_files <- sirius_param(ms1p, result_dir = './MZML/DS_201124_SC_full_PRM_neg_09')


run_sirius <- function(files, ppm_max = 5, ppm_max_ms2 = 15, QC = TRUE, SL = TRUE, SL_path, candidates = 30){
    
    files <- read.csv(files)
    
    for (b in 1:nrow(files)){
        if (QC){
            if (is.na(files[b, "isotopes"])){
                system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                             "formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates,", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"structure --database ALL canopus",
                             sep = " "))
                if(SL){
                    system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNamesSL"],
                             "formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"structure --database",  SL_path ,"canopus",
                             sep = " "))
                }
                
            }
            else if(files[b, "isotopes"] == "present"){
                system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                             "formula --profile orbitrap --candidates", candidates, "--ppm-max", ppm_max,"--ppm-max-ms2", ppm_max_ms2,"structure --database ALL canopus",
                             sep = " "))
                if(SL){
                    system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNamesSL"],
                             "formula --profile orbitrap --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"structure --database",  SL_path ,"canopus",
                             sep = " "))
                }
            }
        }
        else{
            system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                             "formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"structure --database ALL canopus",
                             sep = " "))
            if(SL){
                system(paste("sirius --input", files[b, "sirius_param_file"], "--output", files[b, "outputNamesSL"],
                            "formula --profile orbitrap --no-isotope-filter --no-isotope-score --candidates", candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"structure --database",  SL_path ,"canopus",
                            sep = " "))
            }
        }
    }
}



sirius_postprocess <- function(x, SL = TRUE){
    feat_scale <- function(p) {
    (p - min(p)) / (max(p) - min(p))
    }
    #the result directory name for each file
    dir_name <- paste(x, "/insilico/SIRIUS", sep = '')
    ##result json files in each result directory
    #for suspect list
    SL_Param <- list.files(dir_name, pattern = 'SList.json', full.names = TRUE)
    #for all DB
    Param <- list.files(dir_name, pattern = 'param.json', full.names = TRUE)
    #DATA FRAME of both directories for each feature
    parameter_json <- cbind(SL_Param, Param)
    #read the msdata csv file that contains all features and their data
    msdata <- as.data.frame(read.csv(paste(x, "/insilico/MS1DATA.csv", sep = '')))
    #add new empty columns to store information from SIRIUS results
    msdata$Adducts <- NA           #adducts
    msdata$name <- NA              #name of compound
    msdata$PubChemIDs <- NA        #all pubchem ids
    msdata$SMILES <- NA            #smiles of compound
    msdata$Formula <- NA           #formula of compound
    msdata$FormulaRank <- NA       #formula rank
    msdata$SIRIUSscore <- NA       #Formula score
    msdata$CSIFingerIDscore <- NA  #Structure score
    msdata$SMILESforMCSS <- NA   #SMILES of top scoring candidates; to calculate their tanimoto later 
    msdata$exp_int <- NA           #explained intensity of formula
    msdata$dir <- NA               #directory for results either the strcuture or formula candidates
    msdata$Result <- NA             #the type of candidate, either from Suspect list, strcuture command or formula command
    # for all entries in the json result files for each each feature
    for (i in 1:nrow(parameter_json)){
        
        # find the number of the row corresponding to the row in msdata that has the same precursor m/z as json result files
        rowMS <- msdata[grepl(str_match(as.character(parameter_json[i, 'Param']), "MS1p_\\s*(.*?)\\s*_SIRIUS")[2] ,msdata$premz), ]
        if (SL){
        
            # file path for strcuture candidate from Suspect List json folder
            str_canS <- paste(list.dirs(parameter_json[i,'SL_Param'])[2], '/structure_candidates.tsv', sep = '')
            # file path for formula candidate from Suspect List json folder
            for_canS <- paste(list.dirs(parameter_json[i,'SL_Param'])[2], '/formula_candidates.tsv', sep = '')
            # if the strcuture candidate file exists
            if (file.exists(str_canS)){
            
                # read the corresponding structure and formula candidate files
                str_canSL <- as.data.frame(read_tsv(str_canS))
                for_canSL <- as.data.frame(read_tsv(for_canS))
            
                # if the strcuture candidate file contains 1 or more rows, it has detected a candidate from suspect list, add relevant info
                if (nrow(str_canSL) >= 1){
                
                    # information from structure candidate file
                    msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canSL[1, 'adduct']
                    msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canSL[1, 'name']
                    msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canSL[1, 'pubchemids']
                    msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canSL[1, 'smiles']
                    msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canSL[1, 'molecularFormula']
                    msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canSL[1, 'formulaRank']
                    msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canSL[1, 'CSI:FingerIDScore']
                
                    # information from formula candidate file
                    formulaRow <- which(for_canSL[,'rank'] == str_canSL[1, 'formulaRank'])
                    msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canSL[formulaRow, 'SiriusScore']
                    msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canSL[formulaRow, 'explainedIntensity']
                
                    # other info
                    msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_SL'
                    msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_canS
                }
                # if it's empty, move onto the All DB result folder called PARAM here
                else{
                    # file path for structure and formula candidate from PARAM (ALL DB) json folder
                    str_can <- paste(list.dirs(parameter_json[i,'Param'])[2], '/structure_candidates.tsv', sep = '')
                    for_can <- paste(list.dirs(parameter_json[i,'Param'])[2], '/formula_candidates.tsv', sep = '')
                
                    # if the strcuture candidate file exists
                    if (file.exists(str_can)){
                    
                        # read the corresponding structure and formula candidate files
                        str_canP <- as.data.frame(read_tsv(str_can))
                        for_canP <- as.data.frame(read_tsv(for_can))
                    
                        # if the structure candidate file contains 1 row, it has detected a candidate from all DBs, add relevant info
                        if (nrow(str_canP) == 1){
                        
                            # information from structure candidate file
                            msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP[1, 'adduct']
                            msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP[1, 'name']
                            msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP[1, 'pubchemids']
                            msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP[1, 'smiles']
                            msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP[1, 'molecularFormula']
                            msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP[1, 'formulaRank']
                            msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP[1, 'CSI:FingerIDScore']
                        
                            # information from formula candidate file
                            formulaRow1 <- which(for_canP[,'rank'] == str_canP[1, 'formulaRank'])
                            msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[formulaRow1, 'SiriusScore']
                            msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[formulaRow1, 'explainedIntensity']
                        
                            # other info
                            msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
                            msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can
                        }
                        # if the structure candidate file contains more rows, extract SMILES of top candidates and check their similarity later
                        else if (nrow(str_canP) > 1){
                        
                            # information from structure candidate file
                            msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP[1, 'adduct']
                            msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP[1, 'name']
                            msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP[1, 'pubchemids']
                            msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP[1, 'smiles']
                            msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP[1, 'molecularFormula']
                            msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP[1, 'formulaRank']
                            msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP[1, 'CSI:FingerIDScore']
                        
                            # information from formula candidate file, take info from the formula rank that corresponds to the formula rank with the top strcuture candidate
                            formulaRow2 <- which(for_canP[,'rank'] == str_canP[1, 'formulaRank'])
                            msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[formulaRow2, 'SiriusScore']
                            msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[formulaRow2, 'explainedIntensity']
                        
                            # other info
                            msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
                            msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can
                        
                            # normalize the CSI:FingerIDScores
                            norm_score <- feat_scale(str_canP[,"CSI:FingerIDScore"]) 
                            # store the upper quartile
                            upper_quartile <- str_canP[which(norm_score > as.numeric(quantile(norm_score)[4])), "smiles"]
                        
                            # if the upper quartile has more than 5 candidates, then just take the top 5 candidates
                            if (length(upper_quartile) > 5){
                                upper_quartile <- upper_quartile[1:5]
                            }
                            # save the top candidates SMILES, to check similarity later with rdkit in Python
                            msdata[as.numeric(rownames(rowMS)), 'SMILESforMCSS'] <- paste(upper_quartile, collapse = '|')
                        
                        }
                        # if the structure candidate file is empty, take information from just the formula candidate file
                        else if (nrow(str_canP) == 0){
                        
                            # if formula candidate file is not empty
                            if (nrow(for_canP) >= 1){
                            
                                # information from formula candidate file
                                msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP[1, 'adduct']
                                msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP[1, 'molecularFormula']
                                msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP[1, 'rank']
                                msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[1, 'explainedIntensity']
                                msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[1, 'SiriusScore']
                                # other info
                                msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
                                msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can
                            }
                        }
                    }
                
                    # if the structure candidate from all DBs does not exist
                    else{
                        # check if the formula candidate file exists
                        if (file.exists(for_can)){
                            for_canF1 <- as.data.frame(read_tsv(for_can))
                        
                            # if formula candidate file is not empty
                            if (nrow(for_canF1)>= 1){
                            
                                # information from formula candidate file
                                msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canF1[1, 'adduct']
                                msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canF1[1, 'molecularFormula']
                                msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canF1[1, 'rank']
                                msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canF1[1, 'explainedIntensity']
                                msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canF1[1, 'SiriusScore']
                                # other info
                                msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
                                msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can
                            }
                        }
                    }
                }
            }
            else{
                # directory for structure candidate file 
                str_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/structure_candidates.tsv', sep = '')
                # directory for formula candidate file
                for_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/formula_candidates.tsv', sep = '')
            
                # if the structure candidate file exists
                if (file.exists(str_can1)){
                
                    # read the structure file from All Dbs (PARAM) json file
                    str_canP1 <- as.data.frame(read_tsv(str_can1))
                    # read the formula file from All Dbs (PARAM) json file
                    for_canP1 <- as.data.frame(read_tsv(for_can1))
                
                    #if the structure candidate file has one candidate, it has detected a candidate from all DBs, add relevant info
                    if (nrow(str_canP1) == 1){
                    
                        # information from structure candidate file
                        msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
                        msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
                        msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
                        msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
                        msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
                        msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
                        msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
                        msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
                        # information from formula candidate file
                        formulaRow3 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
                        msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow3, 'SiriusScore']
                        msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow3, 'explainedIntensity']
                        # other info
                        msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
                        msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1
                    
                    }
                    # if the strcuture cabdidate file has more than 1 candidates, it has detected candidates from all DBs, add relevant info
                    else if (nrow(str_canP1) > 1){
                    
                        # information from structure candidate file
                        msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
                        msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
                        msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
                        msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
                        msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
                        msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
                        msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
                        # information from formula candidate file
                        formulaRow4 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
                        msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow4, 'SiriusScore']
                        msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow4, 'explainedIntensity']
                        # other info
                        msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
                        msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1
                    
                        # normalize the CSI:FingerIDScores
                        norm_score1 <- feat_scale(str_canP1[,"CSI:FingerIDScore"]) 
                        # store the upper quartile
                        upper_quartile1 <- str_canP1[which(norm_score1 > as.numeric(quantile(norm_score1)[4])), "smiles"]
                        # if the upper quartile has more than 5 candidates, then just take the top 5 candidates
                        if (length(upper_quartile1) > 5){
                            upper_quartile1 <- upper_quartile1[1:5]
                        }
                        # save the top candidates SMILES, to check similarity later with rdkit in Python
                        msdata[as.numeric(rownames(rowMS)), 'SMILESforMCSS'] <- paste(upper_quartile1, collapse = '|')
                    }
                    # 
                    else if (nrow(str_canP1) == 0){
                        if (file.exists(for_can1)){
                            for_canP2 <- as.data.frame(read_tsv(for_can1))
                            if (nrow(for_canP2)>= 1){
                            
                                # information from formula candidate file
                                msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP2[1, 'adduct']
                                msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP2[1, 'molecularFormula']
                                msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP2[1, 'rank']
                                msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP2[1, 'explainedIntensity']
                                msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP2[1, 'SiriusScore']
                                # other info
                                msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
                                msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
                            }
                        }
                    }
                }
            
                # if the structure candidate file doesn't exists (and no str_sl exists)
                else{
                    # if formula candidate file exists
                    if (file.exists(for_can1)){
                        for_canP3 <- as.data.frame(read_tsv(for_can1))
                        if (nrow(for_canP3) >= 1){
                            # information from formula candidate file
                            msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP3[1, 'adduct']
                            msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP3[1, 'molecularFormula']
                            msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP3[1, 'rank']
                            msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP3[1, 'explainedIntensity']
                            msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP3[1, 'SiriusScore']
                            # other info
                            msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
                            msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
                        }
                    }
                }
            }
        }# SL IS false
        else{
            # directory for structure candidate file 
            str_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/structure_candidates.tsv', sep = '')
            # directory for formula candidate file
            for_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/formula_candidates.tsv', sep = '')
            
            # if the structure candidate file exists
            if (file.exists(str_can1)){
                
                # read the structure file from All Dbs (PARAM) json file
                str_canP1 <- as.data.frame(read_tsv(str_can1))
                # read the formula file from All Dbs (PARAM) json file
                for_canP1 <- as.data.frame(read_tsv(for_can1))
                
                #if the structure candidate file has one candidate, it has detected a candidate from all DBs, add relevant info
                if (nrow(str_canP1) == 1){
                    
                    # information from structure candidate file
                    msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
                    msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
                    msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
                    msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
                    msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
                    msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
                    msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
                    msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
                    # information from formula candidate file
                    formulaRow3 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
                    msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow3, 'SiriusScore']
                    msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow3, 'explainedIntensity']
                    # other info
                    msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
                    msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1
                    
                }
                # if the strcuture cabdidate file has more than 1 candidates, it has detected candidates from all DBs, add relevant info
                else if (nrow(str_canP1) > 1){
                    
                    # information from structure candidate file
                    msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
                    msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
                    msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
                    msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
                    msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
                    msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
                    msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
                    # information from formula candidate file
                    formulaRow4 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
                    msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow4, 'SiriusScore']
                    msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow4, 'explainedIntensity']
                    # other info
                    msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
                    msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1
                    
                    # normalize the CSI:FingerIDScores
                    norm_score1 <- feat_scale(str_canP1[,"CSI:FingerIDScore"]) 
                    # store the upper quartile
                    upper_quartile1 <- str_canP1[which(norm_score1 > as.numeric(quantile(norm_score1)[4])), "smiles"]
                    # if the upper quartile has more than 5 candidates, then just take the top 5 candidates
                    if (length(upper_quartile1) > 5){
                        upper_quartile1 <- upper_quartile1[1:5]
                    }
                    # save the top candidates SMILES, to check similarity later with rdkit in Python
                    msdata[as.numeric(rownames(rowMS)), 'SMILESforMCSS'] <- paste(upper_quartile1, collapse = '|')
                }
                # 
                else if (nrow(str_canP1) == 0){
                    if (file.exists(for_can1)){
                        for_canP2 <- as.data.frame(read_tsv(for_can1))
                        if (nrow(for_canP2)>= 1){
                            
                            # information from formula candidate file
                            msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP2[1, 'adduct']
                            msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP2[1, 'molecularFormula']
                            msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP2[1, 'rank']
                            msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP2[1, 'explainedIntensity']
                            msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP2[1, 'SiriusScore']
                            # other info
                            msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
                            msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
                        }
                    }
                }
            }
            
            # if the structure candidate file doesn't exists (and no str_sl exists)
            else{
                # if formula candidate file exists
                if (file.exists(for_can1)){
                    for_canP3 <- as.data.frame(read_tsv(for_can1))
                    if (nrow(for_canP3) >= 1){
                        # information from formula candidate file
                        msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP3[1, 'adduct']
                        msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP3[1, 'molecularFormula']
                        msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP3[1, 'rank']
                        msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP3[1, 'explainedIntensity']
                        msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP3[1, 'SiriusScore']
                        # other info
                        msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
                        msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
                    }
                }
            }
        }
    }
    write.csv(msdata, paste(x, "/insilico/MS1DATAsirius.csv", sep = ''))
    return(msdata)
}
# Usage: 
# sirius_postprocess(x, SL = TRUE)



#save.image(file = "R_Functions.RData")

#x is the dataframe result from sirius_postprocess
#result_dir 
#input_dir 
#adducts 
#sl_mtfrag 
#SL = TRUE if suspect list present

metfrag_param <- function(x, result_dir, input_dir, adducts, sl_mtfrag, SL = TRUE, ppm_max = 5, ppm_max_ms2= 15){
    
    x <- read.csv(x)
    
    dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/MetFrag", sep =""), "./"), sep = "")
    
    if (!file.exists(dir_name)){
        dir.create(dir_name, recursive = TRUE) ##create folder
    }
    
    db <- c("PubChem", "KEGG")
    
    parameter_file <- c()
    par <- 0
    metfrag_param_file <- c()
    
    AdductsMF <- data.frame(read.csv(adducts))
    
    for (j in 1:nrow(x)){
        if (!(is.na(x[j, 'Adducts']))){
            for (k in db){
                par <- par+1
                para <- as.character(par)
                fileR <- paste(dir_name, "/", para, "_id_", x[j, 'id_X'], "_mz_", x[j, 'premz'], "_rt_", x[j, 'rtmed'], "_db_", k, ".txt", sep = '')
                metfrag_param_file <- c(metfrag_param_file, fileR)

                file.create(fileR, recursive = TRUE)
                file.conn <- file(fileR)
                open(file.conn, open = "at")
                
                
                peakspath <- str_replace(x[j, "ms2Peaks"], "./", input_dir)
                
                
                #writeLines(paste("PeakListPath = ",as.character(peakspath),sep=""),con=file.conn)
                writeLines(paste("PeakListPath = ",peakspath, sep=""),con=file.conn)
                writeLines(paste("IonizedPrecursorMass = ", x[j, "premz"]),con = file.conn)
                
                # write code here
                
                
                PrecursorIonMode <- AdductsMF[which(AdductsMF[, "PrecursorIonType"] == gsub("[[:blank:]]", "", x[j, 'Adducts'])), "PrecursorIonMode"]
                IsPositiveIonMode <- AdductsMF[which(AdductsMF[, "PrecursorIonType"] == gsub("[[:blank:]]", "", x[j, 'Adducts'])), "IsPositiveIonMode"]
            
                
                writeLines(paste("PrecursorIonMode = ", PrecursorIonMode, sep = ''), con = file.conn)
                writeLines(paste("IsPositiveIonMode = ", IsPositiveIonMode, sep = ''), con = file.conn)
                
                writeLines(paste("MetFragDatabaseType = ", k),con = file.conn)
            
                writeLines(paste("DatabaseSearchRelativeMassDeviation = ", ppm_max, sep = ''),con=file.conn)
                writeLines("FragmentPeakMatchAbsoluteMassDeviation = 0.001",con=file.conn)
                writeLines(paste("FragmentPeakMatchRelativeMassDeviation = ", ppm_max_ms2, sep = ''),con=file.conn)
                
                if (SL){
                    writeLines(paste("ScoreSuspectLists = ", sl_mtfrag),con=file.conn)
            
                    writeLines("MetFragScoreTypes = FragmenterScore, SuspectListScore",con=file.conn)
                    writeLines("MetFragScoreWeights = 1.0, 1.0", con=file.conn)
                }
                else{
                    writeLines("MetFragScoreTypes = FragmenterScore", con=file.conn)
                    writeLines("MetFragScoreWeights = 1.0", con=file.conn)
                }
                
                
                writeLines("MetFragCandidateWriter = CSV",con=file.conn)
            
                writeLines(paste("SampleName = ", para, "_id_", x[j, 'id_X'], "_mz_", x[j, 'premz'], "_rt_", x[j, 'rtmed'], "_db_", k, sep = ''),con=file.conn)
                resultspath <- str_replace(result_dir, "./", input_dir)
                writeLines(paste("ResultsPath = ", resultspath, "/insilico/MetFrag/", sep = ''),con=file.conn)
            
                writeLines("MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter",con=file.conn)
                writeLines("MetFragPostProcessingCandidateFilter = InChIKeyFilter",con=file.conn)
                writeLines("MaximumTreeDepth = 2",con=file.conn)
                writeLines("NumberThreads = 1",con=file.conn)
                
                close(file.conn)
                parameter_file <- c(parameter_file,file.conn)
            
            }
        }
    }
    
    write.table(metfrag_param_file, file = paste(result_dir, "/insilico/metparam_list.txt", sep = ""), sep = "/t", row.names = FALSE, col.names = FALSE)
    return(metfrag_param_file)
}

# Usage:
# metfrag_param(x, result_dir, input_dir, adducts, sl_mtfrag, SL = TRUE)



run_metfrag<- function(met_param, input_dir){
    
    filesmet_param <- read.table(met_param)
    
    for (files in filesmet_param[[1]]){
        system(paste("java -jar",  paste(input_dir, "MetFrag2.4.5-CL.jar", sep = ''), files))
    }
}
