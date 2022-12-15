# ---------- Preparations ----------
# Load Libraries
library(Spectra)
library(MsBackendMgf)
library(MsBackendHmdb)
library(MsCoreUtils)
library(MsBackendMsp)
library(readr)
library(dplyr)
# 3 dependencies for latest MassBank version
library(rvest)
library(stringr)
library(xml2)
options(warn=-1)
library("mzR")
#' Download the CompDb database using curl
library(curl)
#library(CompoundDb)

download_specDB_new <- function(input_dir, db = "all"){

    if (dir.exists(input_dir)){
        # Track Time
        start_time <- Sys.time()

        # only input available as of now
        databases <- 'gnps, hmdb, mbank, all'
        # creat a summary file, open and store timings of download and version if possible
        if (!(file.exists(paste(input_dir, "/summaryFile.txt", sep = "")))){
            summaryFile <- paste(input_dir, "/summaryFile.txt", sep = "")
            file.create(summaryFile, recursive = TRUE)
        }
        else{
            summaryFile <- paste(input_dir, "/summaryFile.txt", sep = "")
        }
        file.conn <- file(summaryFile)
        open(file.conn, open = "at")
        # gnps
        if (db == "all" || db =="gnps"){
            # Download file
            system(paste("wget -P", 
                         input_dir, 
                         "https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.msp", 
                         sep =  " "))

            # load the spectra into MsBackendMgf
            gnpsdb <- Spectra(paste(input_dir, "/ALL_GNPS.msp", sep = ''), source = MsBackendMsp())
            save(gnpsdb, file = paste(input_dir,"/gnps.rda", sep = ""))

            # delete the database in its format to free up space
            system(paste("rm", (paste(input_dir, "/ALL_GNPS.mgf", sep = '')), sep = " "))

            writeLines(paste("GNPS saved at", Sys.time(), sep=" "),con=file.conn)
        }
        #mbank
        if (db == "all" || db =="mbank"){

            #print("MassBank WORKS")

            page <- read_html("https://github.com/MassBank/MassBank-data/releases")
            page %>%
                html_nodes("a") %>%       # find all links
                html_attr("href") %>%     # get the url
                str_subset("MassBank_NIST.msp") -> tmp # find those that have the name MassBank_NIST.msp

            #download file
            system(paste("wget ",
                         "https://github.com", tmp[1],
                         sep =  ""))

            mbank <- Spectra(paste(input_dir, "/MassBank_NIST.msp", sep = ''), source = MsBackendMsp())
            save(mbank, file = paste(input_dir,"/mbankNIST.rda", sep = ""))

            # delete the database in its format to free up space
            system(paste("rm", (paste(input_dir, "/MassBank_NIST.msp", sep = '')), sep = " "))

            # obtain the month and year for the database release to add to summary
            res <- str_match(tmp[1], "download/\\s*(.*?)\\s*/MassBank_NIST")

            writeLines(paste("MassBank saved at", Sys.time(), "with release version", res[,2], sep=" "),con=file.conn)
        }

        #mbank
        if (db == "all" || db =="hmdb"){
            # extract HMDB Current version\n",
            html <- read_html("https://hmdb.ca/downloads")
            strings <- html%>% html_elements("a") %>% html_text2()
            ls <- unique(strings)
            hmdb_curr_ver <- c()
            for (i in ls){
                if (grepl("Current", i)){
                hmdb_curr_ver<- c(i, hmdb_curr_ver)
                }
            }
            dbname <- "CompDb.Hsapiens.HMDB.5.0.sqlite"
            db_file <- file.path(tempdir(), dbname)
            curl_download(
                paste0("https://github.com/jorainer/MetaboAnnotationTutorials/",
               "releases/download/2021-11-02/", dbname),
                destfile = db_file)
            #' Load a CompDb database with compound annotation from HMDB
            cdb <- CompDb(db_file)
            hmdb <- Spectra(cdb)
            hmdb$collisionEnergy <- as.numeric(hmdb$collisionEnergy)
            hmdb <- setBackend(hmdb, backend = MsBackendDataFrame())
            save(hmdb, file = paste(input_dir,"/hmdb.rda", sep = ""))
            writeLines(paste("HMDB saved at", Sys.time(), "with release version", hmdb_curr_ver, sep=" "),con=file.conn)
        }

        #wrong input error message
        else if (!grepl(db, databases, fixed = TRUE)){
            stop("Wrong db input. Following inputs apply: gnps, hmdb, mbank or all")
        }
        close(file.conn)
        end_time <- Sys.time()
        print(end_time - start_time)
    }
    else{
        stop("Your input_dir is incorrect. Please provide the directory where all your input files are stored.")
    }
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


##-----------------------------------------------------------------
## normalize intensity 
##-----------------------------------------------------------------

#' Define a function to *normalize* the intensities
norm_int <- function(y, ...) {
    maxint <- max(y[, "intensity"], na.rm = TRUE)
    y[, "intensity"] <- 100 * y[, "intensity"] / maxint
    y
}


## Specifying a function for creating result directories for each input mzml
# input for the function:
# input directory
ms2_rfilename<- function(input_dir){
    if (dir.exists(input_dir)){
        #list_ms2_files <- intersect(list.files(input_dir, pattern = "_PRM_"), list.files(input_dir, pattern = ".mzML"))
        list_ms2_files <- list.files(input_dir, pattern = ".mzML")
        mzml_file <- paste(input_dir, "/", list_ms2_files, sep = "")

        #store the result file names to return to this function as output
        mzml_files <- c()
        ResultFileNames <- c()
        File_id <- c()
        nx <- 0
        # x is mzML files
        for (i in 1:length(mzml_file)){
            nx <- nx+1
            # remove .mzML to extract just the names
            mzml_filex <- str_replace(mzml_file[i], input_dir, ".")
            name_mzmls <- str_remove(as.character(mzml_filex), ".mzML")
            name_mzmlsd <- str_remove(mzml_file[i], ".mzML")
            #name_mzml <- str_replace(name_mzmls, input_dir, "./")
            #' for each file a subdirectory is created to store all results in that, add working directory
            if (!file.exists(name_mzmlsd)){
                dir.create(name_mzmlsd) ##create folder
            }
            ResultFileNames<- c(ResultFileNames, name_mzmls)
            mzml_files <- c(mzml_files, mzml_filex)
            File_id <- c(File_id, paste("file_", nx, sep = ""))
        }
        input_table <- cbind(mzml_files, ResultFileNames, File_id)

        write.csv(input_table, paste(input_dir, "/input_table.csv", sep = ""))
        return(data.frame(input_table))
    }
    else{
        stop("Your input_dir is incorrect. Please provide the directory where all your input files are stored. : ) Good Luck")
    }
}




#' All spectra in mzML files preprocessing, return two outputs, pre-processed MS2 spectra and all precursor masses
# x is one mzML file
spec_Processing <- function(input_dir, x, result_dir){

    x <- paste(input_dir, str_remove(x, "."), sep = "")

    result_dir <- paste(input_dir, str_remove(result_dir, "."), sep = "")
    if (file.exists(x) && substring(x, nchar(x)) == "L"){
        if (dir.exists(result_dir)){
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
            if (!file.exists(paste(result_dir, "/processedSpectra.mzML", sep = ""))){
                export(sps_all, backend = MsBackendMzR(), file = paste(result_dir, "/processedSpectra.mzML", sep = ""))
            }
            if (!file.exists(paste(result_dir, "/premz_list.txt", sep = ""))){
                write.table(pre_mz, file = paste(result_dir, "/premz_list.txt", sep = ""), sep = "/t",row.names = FALSE, col.names = FALSE)
            }

            spsall_pmz <- list(sps_all, pre_mz)
            return(spsall_pmz)
        }
        else{
            stop("Seems like it is not the result directory of the input .mzML file which is provided as x. Please use the function ms2_rfilename to generate a result directory or create one yourself with the same name as the .mzML input file.")
        }
    }
    else{
        stop("Are you sure x is an mzML input file?")
    }

}



spec2_Processing <- function(z, obj, spec = "spec_all", ppmx = 15){
    if (spec == "spec_all"){
        #' Subset the dataset to MS2 spectra matching the m/z
        sps <- filterPrecursorMzValues(obj, mz = z + ppm(c(-z, z), 10))
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
        #' Subset the MB Spectra
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
            m <- m[m[, "mz"] <= z, ]
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


##-----------------------------------------------------------------
## Plotting Mirror Spectra
##-----------------------------------------------------------------

#' Specifying a function to draw peak labels
#label_fun <- function(x) {
    #ints <- unlist(intensity(x))
    #mzs <- format(unlist(mz(x)), digits = 4)
    #mzs[ints < 5] <- ""
    #mzs
#}


spec_dereplication_file <- function(mzml_file, pre_tbl, proc_mzml, db, result_dir, file_id, input_dir, no_of_candidates = 30, ppmx, error = TRUE){
    # if the database selected is HMDB or all
    # if the database selected is GNPS or all
    if (db == "all" || db =="gnps"){
        if (file.exists(paste(input_dir,"/gnps.rda", sep = ""))){
            # load the gnps spectral database
            load(file = paste(input_dir,"/gnps.rda", sep = ""))
        }
        else if (file.exists("gnps.rda")){
            # load the gnps spectral database
            load("gnps.rda")
        }
        
    }
    # if the database selected is HMDB or all
    if (db == "all" || db =="hmdb"){
        
        if (file.exists(paste(input_dir,"/hmdb.rda", sep = ""))){
            # load the hmdb spectral database
            load(file = paste(input_dir,"/hmdb.rda", sep = ""))
        }
        else if (file.exists("hmdb.rda")){
            # load the gnps spectral database
            load("hmdb.rda")
        }
        
    }
    # if the database selected is HMDB or all
    if (db == "all" || db == "mbank"){
        if (file.exists(paste(input_dir,"/mbankNIST.rda", sep = ""))){
            # load the mbank spectral database
            load(file = paste(input_dir,"/mbankNIST.rda", sep = ""))
        }
        else if (file.exists("mbankNIST.rda")){
            # load the gnps spectral database
            load("mbankNIST.rda")
        }
        
    }

    # read spectra object
    sps_all <- Spectra(proc_mzml, source = MsBackendMzR())

    # extract precursor m/z
    tbl <- read.table(pre_tbl)
    pre_mz <- tbl[[1]]

    # common feature information

    id_X <- c() # id
    premz <- c() # precursor mz
    rtmin <- c() # stores rtmin
    rtmax <- c() # stores rtmax
    rtmed <- c() # stores calculated median of rtmin and rtmax
    rtmean <- c() # stores calculated mean of rtmin and rtmax
    col_eng <- c() # stores collision energy
    pol <- c() # stores polarity
    int <- c() # store intensity
    source_file <- c() # source file
    nx <- 0 # numbering the ids

    pre_mzs <- listenv() # list for holding pre_mz futures

    # for each pre mass
    for (x in pre_mz){
        print(x)
        # to name the file
        nx <- nx+1
        # filter spectra based on precusror m/z
        # this is done to extract all common information for id_X
        spsrt <- filterPrecursorMzRange(sps_all, x)
        # id based on file id,
        id_Xx <- paste(file_id,  "M",  as.character(round(x, digits = 0)),
                        "R", as.character(round(median(spsrt$rtime, na.rm = TRUE), digits = 0)),
                        "ID", as.character(nx), sep = '')
        id_X <- c(id_X, id_Xx)
        # pre_mas
        pre <- x
        premz <- c(premz, pre)

        # rt min
        rti <- min(spsrt$rtime)
        rtmin <- c(rtmin, rti)
        #rt max
        rtx <- max(spsrt$rtime)
        rtmax <- c(rtmax, rtx)
        #rt median
        rtmd <- median(spsrt$rtime, na.rm = TRUE)
        rtmed <- c(rtmed, rtmd)
        #rt mean
        rtmn <- mean(spsrt$rtime, na.rm = TRUE)
        rtmean <- c(rtmean, rtmn)
        #collision energy
        ce <- max(spsrt$collisionEnergy)
        col_eng <- c(col_eng, ce)
        #polarity
        pl <- max(spsrt$polarity)
        if (pl == 1){
            px <- 'pos'
            pol <- c(pol, px)
        }
        else {
            px <- 'neg'
            pol <- c(pol, px)
        }
        #int
        ints <- max(spsrt$precursorIntensity)
        int <- c(int, ints)
        #mzmlfile
        source_file <- c(source_file, mzml_file)
        # after all the common infromation is stored,
        # move to extracting matching candidates with input spectra
        pre_mzs[[x]] <- future({
            sps <- spec2_Processing(x, sps_all, spec = "spec_all")
            ####-------------------------------------------------------------
            #### Dereplication with all or GNPS ----
            ####-------------------------------------------------------------
            # define variables for result dataframe

            # if the database selected is GNPS or all
            f_gnps <- future(
            if (db == "all" || db =="gnps"){
                GNPSmax_similarity <- c() # dot product score
                GNPSmzScore <- c() # similar m/z score
                GNPSintScore <- c() # similar int score
                GQMatchingPeaks <- c() # matching peaks between gnps candidate and input spectra
                GNPSTotalPeaks <- c() # total peaks in gnps candidate
                gQueryTotalPeaks<- c() # total peaks in input spectra
                GNPSSMILES <- c() # smiles of gnps candidate
                #GNPSspectrumID <- c() # spectrum id of gnps candidate
                GNPScompound_name <- c() # compound name of gnps candidate
                #GNPSmirrorSpec <- c() # path for mirror spectra between gnps candidate and input of gnps candidate
                Source <- c() # GNPS as source of result
                #### GNPS spec with pre_mz
                gnps_with_mz <- spec2_Processing(x, gnpsdb, spec = "gnps", ppmx) # change here later

                # define the directoyr name to store all GNPS results
                dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/GNPS/", sep = ""), "."), sep ="")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
                if (length(sps) != 0 && length(gnps_with_mz) !=0){
                     #' Compare experimental spectra against GNPS
                    res <- compareSpectra(sps, gnps_with_mz, ppm = 15, FUN = MsCoreUtils::gnps, MAPFUN = joinPeaksGnps)

                    # first condition for GNPS
                    # if more input spectra and more candidates have been extracted from GNPS
                    if (length(sps) > 1 && length(gnps_with_mz) >1){
                        # given threshold of 0.85 for GNPS, extract top candidates
                        res_top <- which(res > res[res>0.85], arr.ind = TRUE)

                        # if there are some compounds from GNPS detected
                        if (length(res_top) > 0){
                            res_topdf <- data.frame(res_top)
                            # to store the scores to add to res_topdf
                            gnps_scores <- c()

                            # for all rows and columns in res_topdf
                            for (i in 1:nrow(res_topdf)){

                                # store the scores
                                gnps_scores <- c(gnps_scores, res[(res_topdf[i, "row"]), (res_topdf[i, "col"])])
                            }
                            if (length(gnps_scores)>0){
                                # add the score column to res_top
                                gnps_res <- cbind(res_top, gnps_scores)
                                gnps_res <- data.frame(gnps_res)

                                # sort in descending order
                                ordered_gnps_res <- gnps_res[order(-gnps_res[,"gnps_scores"]),]
                                df_ord_gnps_res <- data.frame(ordered_gnps_res)
                                if (nrow(df_ord_gnps_res)>no_of_candidates){
                                    df_ord_gnps_res <- df_ord_gnps_res[1:no_of_candidates,]
                                }
                                #for each candidate from GNPS
                                for (k in 1:nrow(df_ord_gnps_res)){
                                    # take each component from df_ord_gnps_res
                                    idv <- df_ord_gnps_res[k,]
                                    df_peaklists <- peakdf(gnps_with_mz[idv[[2]]], sps[idv[[1]]], ppmx)

                                    if (!(is.null(df_peaklists))){
                                        GNPSscore <- idv[1, "gnps_scores"]
                                        GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)

                                        GNPSmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_with_mz[idv[[2]]])[[1]])+nrow(peaksData(sps[idv[[1]]])[[1]]))
                                        GNPSmzScore <- c(GNPSmzScore, GNPSmz)

                                        GNPSint <- mean(1-(df_peaklists[,"diff"]/100))
                                        GNPSintScore <- c(GNPSintScore, GNPSint)

                                        GQMatPeaks <- nrow(df_peaklists)
                                        GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)

                                        GNPSTPeaks <- nrow(peaksData(gnps_with_mz[idv[[2]]])[[1]])
                                        GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)

                                        gQTPeaks<- nrow(peaksData(sps[idv[[1]]])[[1]])
                                        gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)

                                        GNPS_SMILES <- gnps_with_mz[idv[[2]]]$smiles
                                        GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)

                                        #GNPSID <- gnps_with_mz[idv[[2]]]$SPECTRUMID
                                        #GNPSspectrumID <- c(GNPSspectrumID, GNPSID)

                                        GNPSname <- gnps_with_mz[idv[[2]]]$name
                                        GNPScompound_name <- c(GNPScompound_name, GNPSname)
                                        Src <- "GNPS"
                                        Source <- c(Source, Src)

                                    }# if df_peaklists isnt empty
                                }# for each candidate

                            }# gnps_score exists
                        }# if res_top has some good candidates
                    }# first condition
                    # if only one sepctrum from input and more candidates from GNPS
                    else if (length(sps) == 1 && length(gnps_with_mz) >1){
                        # given threshold of 0.85 for GNPS, extract top candidates
                        res_top <- which(res > res[res>0.85], arr.ind = TRUE)
                        # if there are candidates with good score
                        if (length(res_top) > 0){
                            res_topdf <- data.frame(res_top)

                            # top store the scores to add to res_topdf
                            gnps_scores <- c()
                            # for all rows and columns in res_topdf

                            for (i in 1:nrow(res_topdf)){
                                # store the scores
                                gnps_scores <- c(gnps_scores, res[(res_topdf[i, "res_topdf"])])
                            }
                            if (length(gnps_scores)>0){
                                # add the score column to res_top
                                gnps_res <- cbind(res_top, gnps_scores)
                                gnps_res <- data.frame(gnps_res)

                                # sort in descending order
                                ordered_gnps_res <- gnps_res[order(-gnps_res[,"gnps_scores"]),]

                                df_ord_gnps_res <- data.frame(ordered_gnps_res)
                                if (nrow(df_ord_gnps_res)>no_of_candidates){
                                    df_ord_gnps_res <- df_ord_gnps_res[1:no_of_candidates,]
                                }
                                # for each candidate
                                for (k in 1:nrow(df_ord_gnps_res)){
                                    # take each candidate
                                    idv <- df_ord_gnps_res[k,]
                                    df_peaklists <- peakdf(gnps_with_mz[idv[[1]]], sps, ppmx)

                                    # if there are matchingpeaks
                                    if (!(is.null(df_peaklists))){
                                        GNPSscore <- idv[1, "gnps_scores"]
                                        GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)

                                        GNPSmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_with_mz[idv[[1]]])[[1]])+nrow(peaksData(sps)[[1]]))
                                        GNPSmzScore <- c(GNPSmzScore, GNPSmz)

                                        GNPSint <- mean(1-(df_peaklists[,"diff"]/100))
                                        GNPSintScore <- c(GNPSintScore, GNPSint)

                                        GQMatPeaks <- nrow(df_peaklists)
                                        GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)

                                        GNPSTPeaks <- nrow(peaksData(gnps_with_mz[idv[[1]]])[[1]])
                                        GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)

                                        gQTPeaks<- nrow(peaksData(sps)[[1]])
                                        gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)

                                        GNPS_SMILES <- gnps_with_mz[idv[[1]]]$smiles
                                        GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)

                                        #GNPSID <- gnps_with_mz[idv[[1]]]$SPECTRUMID
                                        #GNPSspectrumID <- c(GNPSspectrumID, GNPSID)

                                        GNPSname <- gnps_with_mz[idv[[1]]]$name
                                        GNPScompound_name <- c(GNPScompound_name, GNPSname)

                                        Src <- "GNPS"
                                        Source <- c(Source, Src)
                                    }# if df_peaklists isnt empty
                                }# for each candidate
                            }# gnps_score exists
                        }# if res_top has some good candidates
                    }# second condition
                    # if there are more input spectra and one candidate from GNPS
                    else if (length(sps) > 1 && length(gnps_with_mz) == 1){
                        # given threshold of 0.85 for GNPS, extract top candidates
                        res_top <- which(res > res[res>0.85], arr.ind = TRUE)

                        # if there are good matching candidates
                        if (length(res_top) > 0){
                            res_topdf <- data.frame(res_top)

                            # top store the scores to add to res_topdf
                            gnps_scores <- c()
                            # for all rows and columns in res_topdf

                            # for all candidates
                            for (i in 1:nrow(res_topdf)){
                                # store the scores
                                gnps_scores <- c(gnps_scores, res[(res_topdf[i, "res_topdf"])])
                            }
                            if (length(gnps_scores)>0){
                                # add the score column to res_top
                                gnps_res <- cbind(res_top, gnps_scores)
                                gnps_res <- data.frame(gnps_res)

                                # sort in descending order
                                ordered_gnps_res <- gnps_res[order(-gnps_res[,"gnps_scores"]),]

                                df_ord_gnps_res <- data.frame(ordered_gnps_res)
                                if (nrow(df_ord_gnps_res)>no_of_candidates){
                                    df_ord_gnps_res <- df_ord_gnps_res[1:no_of_candidates,]
                                }
                                # for each candidate match
                                for (k in 1:nrow(df_ord_gnps_res)){
                                    # take each candidate
                                    idv <- df_ord_gnps_res[k,]

                                    df_peaklists <- peakdf(gnps_with_mz, sps[idv[[1]]], ppmx)

                                    # if there are matching peaks
                                    if (!(is.null(df_peaklists))){
                                        GNPSscore <- idv[1, "gnps_scores"]
                                        GNPSmax_similarity <- c(GNPSmax_similarity, GNPSscore)

                                        GNPSmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(gnps_with_mz)[[1]])+nrow(peaksData(sps[idv[[1]]])[[1]]))
                                        GNPSmzScore <- c(GNPSmzScore, GNPSmz)

                                        GNPSint <- mean(1-(df_peaklists[,"diff"]/100))
                                        GNPSintScore <- c(GNPSintScore, GNPSint)

                                        GQMatPeaks <- nrow(df_peaklists)
                                        GQMatchingPeaks <- c(GQMatchingPeaks, GQMatPeaks)

                                        GNPSTPeaks <- nrow(peaksData(gnps_with_mz)[[1]])
                                        GNPSTotalPeaks <- c(GNPSTotalPeaks, GNPSTPeaks)

                                        gQTPeaks<- nrow(peaksData(sps[idv[[1]]])[[1]])
                                        gQueryTotalPeaks <- c(gQueryTotalPeaks, gQTPeaks)

                                        GNPS_SMILES <- gnps_with_mz$smiles
                                        GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)

                                        #GNPSID <- gnps_with_mz$SPECTRUMID
                                        #GNPSspectrumID <- c(GNPSspectrumID, GNPSID)

                                        GNPSname <- gnps_with_mz$name
                                        GNPScompound_name <- c(GNPScompound_name, GNPSname)

                                        Src <- "GNPS"
                                        Source <- c(Source, Src)
                                    }# if df_peaklists isnt empty
                                }# for each candidate
                            }# gnps_score exists 
                        }# if res_top has some good candidates
                    }# third condition
                    else if (length(sps) == 1 && length(gnps_with_mz) == 1){
                        if (res>= 0.85){
                            #take that one candidate
                            gnps_best_match <- gnps_with_mz

                            df_peaklists <- peakdf(gnps_best_match, sps, ppmx)

                            # if there are matching peaks
                            if (!(is.null(df_peaklists))){

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


                                GNPS_SMILES <- gnps_best_match$smiles
                                GNPSSMILES <- c(GNPSSMILES, GNPS_SMILES)

                                GNPSname <- gnps_best_match$name
                                GNPScompound_name <- c(GNPScompound_name, GNPSname)


                                #GNPSID <- gnps_best_match$SPECTRUMID
                                #GNPSspectrumID <- c(GNPSspectrumID, GNPSID)


                                Src <- "GNPS"
                                Source <- c(Source, Src)
                            }# if candidate exixts
                        }#  if res_top has some good candidates

                    }# fourth condition
                }# if sps and gnps has some candidate
                gnps_x <- data.frame(cbind(GNPSmax_similarity, GNPSmzScore,
                                       GNPSintScore, GQMatchingPeaks,
                                       GNPSTotalPeaks, gQueryTotalPeaks,
                                       GNPSSMILES, GNPScompound_name, Source))
                write.csv(gnps_x, file = paste(dir_name, "/gnps_results_for_", id_Xx, ".csv", sep = ""))
            })# gnps ends


            ####-------------------------------------------------------------
            #### Dereplication with all or HMDB ----
            ####-------------------------------------------------------------
             # if the database selected is HMDB or all
            f_hmdb <- future(
            if (db == "all" || db =="hmdb"){
                # hmdb
                HMDBmax_similarity <- c()
                HMDBmzScore <- c()
                HMDBintScore <- c()
                HQMatchingPeaks <- c()
                HMDBTotalPeaks <- c()
                hQueryTotalPeaks<- c()
                HMDBcompoundID <- c()
                Source <- c()

                #### HMDB spec with pre_mz
                hmdb_with_mz <- spec2_Processing(x, hmdb, spec = "hmdb", ppmx) # change here later

                # directory name for HMDB results
                dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/HMDB/", sep = ""), "."), sep ="")
                if (!file.exists(dir_name)){
                    dir.create(dir_name, recursive = TRUE)
                }
                if (length(sps) != 0 && length(hmdb_with_mz) !=0){
                    #' Compare experimental spectra against HMDB
                    res <- compareSpectra(sps, hmdb_with_mz, ppm = 15)

                    # if there are more input spectra and more candidates from GNPS
                    if (length(sps) > 1 && length(hmdb_with_mz) >1){
                        # given threshold of 0.70 for HMDB, extract top candidates
                        res_top <- which(res > res[res>0.70], arr.ind = TRUE)
                        if (length(res_top) > 0){
                            res_topdf <- data.frame(res_top)
                            # to store the scores to add to res_topdf
                            hmdb_scores <- c()
                            # for all rows and columns in res_topdf
                            for (i in 1:nrow(res_topdf)){

                                # store the scores
                                hmdb_scores <- c(hmdb_scores, res[(res_topdf[i, "row"]), (res_topdf[i, "col"])])
                            }
                            if (length(hmdb_scores)>0){
                                # add the score column to res_top
                                hmdb_res <- cbind(res_top, hmdb_scores)
                                hmdb_res <- data.frame(hmdb_res)
                                # sort in descending order
                                ordered_hmdb_res <- hmdb_res[order(-hmdb_res[,"hmdb_scores"]),]
                                df_ord_hmdb_res <- data.frame(ordered_hmdb_res)
                                if (nrow(df_ord_hmdb_res)>no_of_candidates){
                                    df_ord_hmdb_res <- df_ord_hmdb_res[1:no_of_candidates,]
                                }
                                for (k in 1:nrow(df_ord_hmdb_res)){
                                    idv <- df_ord_hmdb_res[k,]
                                    df_peaklists <- peakdf(hmdb_with_mz[idv[[2]]], sps[idv[[1]]], ppmx)
                                    if (!(is.null(df_peaklists))){
                                        HMDBscore <- idv[1, "hmdb_scores"]
                                        HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)

                                        HMDBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_with_mz[idv[[2]]])[[1]])+nrow(peaksData(sps[idv[[1]]])[[1]]))
                                        HMDBmzScore <- c(HMDBmzScore, HMDBmz)

                                        HMDBint <- mean(1-(df_peaklists[,"diff"]/100))
                                        HMDBintScore <- c(HMDBintScore, HMDBint)


                                        HQMatPeaks <- nrow(df_peaklists)
                                        HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                                        HMDBTPeaks <- nrow(peaksData(hmdb_with_mz[idv[[2]]])[[1]])
                                        HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)

                                        hQTPeaks<- nrow(peaksData(sps[idv[[1]]])[[1]])
                                        hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                                        HMDBID <- hmdb_with_mz[idv[[2]]]$compound_id
                                        HMDBcompoundID <- c(HMDBcompoundID, HMDBID)

                                        Src <- "HMDB"
                                        Source <- c(Source, Src)
                                    }# if df_peaklists is not empty
                                }# for each candidate
                            }# if hmdb_scores exist
                        }# if there are top candidadates with good scores
                    }#first condition ends
                    else if (length(sps) == 1 && length(hmdb_with_mz) >1){
                        # given threshold of 0.70 for HMDB, extract top candidates
                        res_top <- which(res > res[res>0.70], arr.ind = TRUE)
                        if (length(res_top) > 0){
                            res_topdf <- data.frame(res_top)

                            # top store the scores to add to res_topdf
                            hmdb_scores <- c()
                            # for all rows and columns in res_topdf
                            for (i in 1:nrow(res_topdf)){
                                # store the scores
                                hmdb_scores <- c(hmdb_scores, res[(res_topdf[i, "res_topdf"])])
                            }
                            if (length(hmdb_scores > 0)){
                                # add the score column to res_top
                                hmdb_res <- cbind(res_top, hmdb_scores)
                                hmdb_res <- data.frame(hmdb_res)
                                # sort in descending order
                                ordered_hmdb_res <- hmdb_res[order(-hmdb_res[,"hmdb_scores"]),]

                                df_ord_hmdb_res <- data.frame(ordered_hmdb_res)
                                if (nrow(df_ord_hmdb_res)>no_of_candidates){
                                    df_ord_hmdb_res <- df_ord_hmdb_res[1:no_of_candidates,]
                                }
                                for (k in 1:nrow(df_ord_hmdb_res)){
                                    idv <- df_ord_hmdb_res[k,]
                                    df_peaklists <- peakdf(hmdb_with_mz[idv[[1]]], sps, ppmx)
                                    if (!(is.null(df_peaklists))){
                                        HMDBscore <- idv[1, "hmdb_scores"]
                                        HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)

                                        HMDBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_with_mz[idv[[1]]])[[1]])+nrow(peaksData(sps)[[1]]))
                                        HMDBmzScore <- c(HMDBmzScore, HMDBmz)

                                        HMDBint <- mean(1-(df_peaklists[,"diff"]/100))
                                        HMDBintScore <- c(HMDBintScore, HMDBint)


                                        HQMatPeaks <- nrow(df_peaklists)
                                        HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                                        HMDBTPeaks <- nrow(peaksData(hmdb_with_mz[idv[[1]]])[[1]])
                                        HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)

                                        hQTPeaks<- nrow(peaksData(sps)[[1]])
                                        hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                                        HMDBID <- hmdb_with_mz[idv[[1]]]$compound_id
                                        HMDBcompoundID <- c(HMDBcompoundID, HMDBID)

                                        Src <- "HMDB"
                                        Source <- c(Source, Src)
                                    }#df_peaklists is not null
                                }# for each candidate
                            }# hmdb_scores exists
                        }
                    }# second condition ends
                    else if (length(sps) > 1 && length(hmdb_with_mz) == 1){
                        # given threshold of 0.70 for HMDB, extract top candidates
                        res_top <- which(res > res[res>0.70], arr.ind = TRUE)
                        if (length(res_top) > 0){
                            res_topdf <- data.frame(res_top)
                            # top store the scores to add to res_topdf
                            hmdb_scores <- c()
                            # for all rows and columns in res_topdf

                            for (i in 1:nrow(res_topdf)){
                                # store the scores
                                hmdb_scores <- c(hmdb_scores, res[(res_topdf[i, "res_topdf"])])
                            }
                            if (length(hmdb_scores)>0){
                                # add the score column to res_top
                                hmdb_res <- cbind(res_top, hmdb_scores)
                                hmdb_res <- data.frame(hmdb_res)

                                # sort in descending order
                                ordered_hmdb_res <- hmdb_res[order(-hmdb_res[,"hmdb_scores"]),]

                                df_ord_hmdb_res <- data.frame(ordered_hmdb_res)
                                if (nrow(df_ord_hmdb_res)>no_of_candidates){
                                    df_ord_hmdb_res <- df_ord_hmdb_res[1:no_of_candidates,]
                                }
                                for (k in 1:nrow(df_ord_hmdb_res)){
                                    idv <- df_ord_hmdb_res[k,]
                                    df_peaklists <- peakdf(hmdb_with_mz, sps[idv[[1]]], ppmx)
                                    if (!(is.null(df_peaklists))){
                                        HMDBscore <- idv[1, "hmdb_scores"]
                                        HMDBmax_similarity <- c(HMDBmax_similarity, HMDBscore)

                                        HMDBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(hmdb_with_mz)[[1]])+nrow(peaksData(sps[idv[[1]]])[[1]]))
                                        HMDBmzScore <- c(HMDBmzScore, HMDBmz)

                                        HMDBint <- mean(1-(df_peaklists[,"diff"]/100))
                                        HMDBintScore <- c(HMDBintScore, HMDBint)


                                        HQMatPeaks <- nrow(df_peaklists)
                                        HQMatchingPeaks <- c(HQMatchingPeaks, HQMatPeaks)


                                        HMDBTPeaks <- nrow(peaksData(hmdb_with_mz)[[1]])
                                        HMDBTotalPeaks <- c(HMDBTotalPeaks, HMDBTPeaks)

                                        hQTPeaks<- nrow(peaksData(sps[idv[[1]]])[[1]])
                                        hQueryTotalPeaks<- c(hQueryTotalPeaks, hQTPeaks)


                                        HMDBID <- hmdb_with_mz$compound_id
                                        HMDBcompoundID <- c(HMDBcompoundID, HMDBID)

                                        Src <- "HMDB"
                                        Source <- c(Source, Src)
                                    }# if df_peaklists isn't null
                                }# for each candidate
                            }#hmdb_score exists
                        }
                    }#third condition ends
                    else if (length(sps) == 1 && length(hmdb_with_mz) == 1){
                        if (res>=0.70){
                            hmdb_best_match <- hmdb_with_mz
                            df_peaklists <- peakdf(hmdb_best_match, sps, ppmx)
                            if (!(is.null(df_peaklists))){
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
                                 Src <- "HMDB"
                                Source <- c(Source, Src)
                            }
                        }
                    }# fourth condition
                }# if both sps and hmdb has some matching candidates
                hmdb_x <- data.frame(cbind(HMDBmax_similarity, HMDBmzScore,
                                       HMDBintScore, HQMatchingPeaks,
                                       HMDBTotalPeaks, hQueryTotalPeaks,
                                       HMDBcompoundID, Source))
                write.csv(hmdb_x, file = paste(dir_name, "/hmdb_results_for_", id_Xx, ".csv", sep = ""))
            })#hmdb ends here



            ####-------------------------------------------------------------
            #### Dereplication with all or MassBank ----
            ####-------------------------------------------------------------
            # define variables for result dataframe
            f_mbank <- future(
                # if the database selected is MassBank or all
                if (db == "all" || db =="mbank"){
                    # mbank
                    MBmax_similarity <- c()
                    MBmzScore <- c()
                    MBintScore <- c()
                    MQMatchingPeaks <- c()
                    MBTotalPeaks <- c()
                    mQueryTotalPeaks<- c()
                    MBformula <- c()
                    MBSMILES <- c()
                    MBspectrumID <- c()
                    MBcompound_name <- c()
                    Source <- c()

                    mbank_with_mz <- spec2_Processing(x, mbank, spec = "mbank", ppmx) # change here later

                    dir_name <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication/MassBank/", sep = ""), "."), sep ="")
                    if (!file.exists(dir_name)){
                        dir.create(dir_name, recursive = TRUE)
                    }
                    if (length(sps) != 0 && length(mbank_with_mz) !=0){
                        #' Compare experimental spectra against MassBank
                        res <- compareSpectra(sps, mbank_with_mz, ppm = 15)
                        if (length(sps) > 1 && length(mbank_with_mz) >1){
                            # given threshold of 0.70 for MassBank, extract top candidates
                            res_top <- which(res > res[res>0.70], arr.ind = TRUE)
                            if (length(res_top) > 0){
                                res_topdf <- data.frame(res_top)
                                # to store the scores to add to res_topdf
                                mbank_scores <- c()
                                # for all rows and columns in res_topdf
                                for (i in 1:nrow(res_topdf)){
                                    # store the scores
                                    mbank_scores <- c(mbank_scores, res[(res_topdf[i, "row"]), (res_topdf[i, "col"])])
                                }
                                if (length(mbank_scores)>0){
                                    # add the score column to res_top
                                    mbank_res <- cbind(res_top, mbank_scores)
                                    mbank_res <- data.frame(mbank_res)
                                    # sort in descending order
                                    ordered_mbank_res <- mbank_res[order(-mbank_res[,"mbank_scores"]),]
                                    df_ord_mbank_res <- data.frame(ordered_mbank_res)
                                    if (nrow(df_ord_mbank_res)>no_of_candidates){
                                        df_ord_mbank_res <- df_ord_mbank_res[1:no_of_candidates,]
                                    }
                                    for (k in 1:nrow(df_ord_mbank_res)){
                                        idv <- df_ord_mbank_res[k,]
                                        df_peaklists <- peakdf(mbank_with_mz[idv[[2]]], sps[idv[[1]]], ppmx)
                                        if (!(is.null(df_peaklists))){
                                            mbscore <- idv[1, "mbank_scores"]
                                            MBmax_similarity<- c(MBmax_similarity, mbscore)

                                            MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_with_mz[idv[[2]]])[[1]])+nrow(peaksData(sps[idv[[1]]])[[1]]))
                                            MBmzScore <- c(MBmzScore, MBmz)

                                            MBint <- mean(1-(df_peaklists[,"diff"]/100))
                                            MBintScore <- c(MBintScore, MBint)

                                            MQMatPeaks <- nrow(df_peaklists)
                                            MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                                            MBTPeaks <- nrow(peaksData(mbank_with_mz[idv[[2]]])[[1]])
                                            MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                                            mQTPeaks<- nrow(peaksData(sps[idv[[1]]])[[1]])
                                            mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)


                                            MBfor <- mbank_with_mz[idv[[2]]]$Formula
                                            MBformula<- c(MBformula, MBfor)

                                            MBS <- mbank_with_mz[idv[[2]]]$smiles
                                            MBSMILES <- c(MBSMILES, MBS)

                                            MBID <- mbank_with_mz[idv[[2]]]$accession
                                            MBspectrumID<- c(MBspectrumID, MBID)

                                            MBname <- mbank_with_mz[idv[[2]]]$Name
                                            MBcompound_name <- c(MBcompound_name, MBname)

                                            Src <- "MassBank"
                                            Source <- c(Source, Src)
                                        }# if df_peaklists is not empty
                                    }# for each candidate
                                }# if mbank_score exists
                            }# if there are in candidates in res_top
                        }# first condition
                        else if (length(sps) == 1 && length(mbank_with_mz) >1){
                            # given threshold of 0.70 for MassBank, extract top candidates
                            res_top <- which(res > res[res>0.70], arr.ind = TRUE)
                            if (length(res_top) > 0){
                                res_topdf <- data.frame(res_top)

                                # top store the scores to add to res_topdf
                                mbank_scores <- c()
                                # for all rows and columns in res_topdf

                                for (i in 1:nrow(res_topdf)){
                                    # store the scores
                                    mbank_scores <- c(mbank_scores, res[(res_topdf[i, "res_top"])])
                                }
                                if (length(mbank_scores)>0){
                                    # add the score column to res_top
                                    mbank_res <- cbind(res_top, mbank_scores)
                                    mbank_res <- data.frame(mbank_res)

                                    # sort in descending order
                                    ordered_mbank_res <- mbank_res[order(-mbank_res[,"mbank_scores"]),]

                                    df_ord_mbank_res <- data.frame(ordered_mbank_res)
                                    if (nrow(df_ord_mbank_res)>no_of_candidates){
                                        df_ord_mbank_res <- df_ord_mbank_res[1:no_of_candidates,]
                                    }

                                    for (k in 1:nrow(df_ord_mbank_res)){
                                        idv <- df_ord_mbank_res[k,]
                                        df_peaklists <- peakdf(mbank_with_mz[idv[[1]]], sps, ppmx)
                                        if (!(is.null(df_peaklists))){
                                            mbscore <- idv[1, "mbank_scores"]
                                            MBmax_similarity<- c(MBmax_similarity, mbscore)

                                            MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_with_mz[idv[[1]]])[[1]])+nrow(peaksData(sps)[[1]]))
                                            MBmzScore <- c(MBmzScore, MBmz)

                                            MBint <- mean(1-(df_peaklists[,"diff"]/100))
                                            MBintScore <- c(MBintScore, MBint)

                                            MQMatPeaks <- nrow(df_peaklists)
                                            MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                                            MBTPeaks <- nrow(peaksData(mbank_with_mz[idv[[1]]])[[1]])
                                            MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                                            mQTPeaks<- nrow(peaksData(sps)[[1]])
                                            mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)


                                            MBfor <- mbank_with_mz[idv[[1]]]$Formula
                                            MBformula<- c(MBformula, MBfor)

                                            MBS <-  mbank_with_mz[idv[[1]]]$smiles
                                            MBSMILES <- c(MBSMILES, MBS)

                                            MBID <- mbank_with_mz[idv[[1]]]$accession
                                            MBspectrumID<- c(MBspectrumID, MBID)

                                            MBname <- mbank_with_mz[idv[[1]]]$Name
                                            MBcompound_name <- c(MBcompound_name, MBname)

                                            Src <- "MassBank"
                                            Source <- c(Source, Src)
                                        }# if df_peaklists is not empty
                                    }# for each candidate
                                }# if mbank_score exists

                            }# res_top has some candidates
                        }#second condition ends
                        else if (length(sps) > 1 && length(mbank_with_mz) == 1){
                            # given threshold of 0.70 for MassBank, extract top candidates
                            res_top <- which(res > res[res>0.70], arr.ind = TRUE)
                            if (length(res_top) > 0){
                                res_topdf <- data.frame(res_top)

                                # top store the scores to add to res_topdf
                                mbank_scores <- c()
                                # for all rows and columns in res_topdf

                                for (i in 1:nrow(res_topdf)){
                                    # store the scores
                                    mbank_scores <- c(mbank_scores, res[(res_topdf[i, "res_top"])])
                                }
                                if (length(mbank_scores)>0){
                                    # add the score column to res_top
                                    mbank_res <- cbind(res_top, mbank_scores)
                                    mbank_res <- data.frame(mbank_res)
                                    # sort in descending order
                                    ordered_mbank_res <- mbank_res[order(-mbank_res[,"mbank_scores"]),]

                                    df_ord_mbank_res <- data.frame(ordered_mbank_res)
                                    if (nrow(df_ord_mbank_res)>no_of_candidates){
                                        df_ord_mbank_res <- df_ord_mbank_res[1:no_of_candidates,]
                                    }
                                    for (k in 1:nrow(df_ord_mbank_res)){
                                        idv <- df_ord_mbank_res[k,]
                                        df_peaklists <- peakdf(mbank_with_mz, sps[idv[[1]]], ppmx)
                                        if (!(is.null(df_peaklists))){
                                            mbscore <- idv[1, "mbank_scores"]
                                            MBmax_similarity<- c(MBmax_similarity, mbscore)

                                            MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_with_mz)[[1]])+nrow(peaksData(sps[idv[[1]]])[[1]]))
                                            MBmzScore <- c(MBmzScore, MBmz)

                                            MBint <- mean(1-(df_peaklists[,"diff"]/100))
                                            MBintScore <- c(MBintScore, MBint)

                                            MQMatPeaks <- nrow(df_peaklists)
                                            MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                                            MBTPeaks <- nrow(peaksData(mbank_with_mz)[[1]])
                                            MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                                            mQTPeaks<- nrow(peaksData(sps[idv[[1]]])[[1]])
                                            mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)


                                            MBfor <- mbank_with_mz$Formula
                                            MBformula<- c(MBformula, MBfor)

                                            MBS <- mbank_with_mz$smiles
                                            MBSMILES <- c(MBSMILES, MBS)
                                            MBID <- mbank_with_mz$accession
                                            MBspectrumID<- c(MBspectrumID, MBID)

                                            MBname <- mbank_with_mz$Name
                                            MBcompound_name <- c(MBcompound_name, MBname)

                                            Src <- "MassBank"
                                            Source <- c(Source, Src)
                                        }# if df_peaklists is not empty
                                    }# for each candidate
                                }# if mbank_score exists
                            }# res_top has some candidates
                        }#third condition ends
                        else if (length(sps) == 1 && length(mbank_with_mz) == 1){
                            if (res>0.70){
                                mbank_best_match <- mbank_with_mz
                                df_peaklists <- peakdf(mbank_best_match, sps, ppmx)
                                if (!(is.null(df_peaklists))){
                                    mbscore <- max(res)
                                    MBmax_similarity<- c(MBmax_similarity, mbscore)

                                    MBmz <- (nrow(df_peaklists)*2)/(nrow(peaksData(mbank_with_mz)[[1]])+nrow(peaksData(sps)[[1]]))
                                    MBmzScore <- c(MBmzScore, MBmz)

                                    MBint <- mean(1-(df_peaklists[,"diff"]/100))
                                    MBintScore <- c(MBintScore, MBint)

                                    MQMatPeaks <- nrow(df_peaklists)
                                    MQMatchingPeaks <- c(MQMatchingPeaks, MQMatPeaks)

                                    MBTPeaks <- nrow(peaksData(mbank_with_mz)[[1]])
                                    MBTotalPeaks<- c(MBTotalPeaks, MBTPeaks)

                                    mQTPeaks<- nrow(peaksData(sps)[[1]])
                                    mQueryTotalPeaks<- c(mQueryTotalPeaks, mQTPeaks)


                                    MBfor <- mbank_with_mz$Formula
                                    MBformula<- c(MBformula, MBfor)

                                    MBS <- mbank_with_mz$smiles
                                    MBSMILES <- c(MBSMILES, MBS)

                                    MBID <- mbank_with_mz$accession
                                    MBspectrumID<- c(MBspectrumID, MBID)

                                    MBname <- mbank_with_mz$Name
                                    MBcompound_name <- c(MBcompound_name, MBname)

                                    Src <- "MassBank"
                                    Source <- c(Source, Src)

                                }# if df_peaklists is not empty
                            }# res_top has some candidates
                        }#fourth condition ends
                    }# both sps and mbank have candidates
                    mbank_x <- data.frame(cbind(MBmax_similarity, MBmzScore,
                                            MBintScore, MQMatchingPeaks,
                                            MBTotalPeaks, mQueryTotalPeaks,
                                            MBformula, MBSMILES, MBspectrumID,
                                            MBcompound_name, Source))
                    write.csv(mbank_x, file = paste(dir_name, "/mbank_results_for_", id_Xx, ".csv", sep = ""))

                }
            )
            v <- c(future::value(f_gnps), future::value(f_hmdb), future::value(f_mbank))

        })
    }
    pre_mzs <- as.list(pre_mzs)
    v_pre_mzs <- future::value(pre_mzs)
    result_dir_spectra <- paste(input_dir, str_remove(paste(result_dir, "/spectral_dereplication", sep = ""), "."), sep = "")
    spectra_input <- data.frame(cbind(id_X, premz, rtmin,
                                      rtmax, rtmed, rtmean,
                                      col_eng, pol, int, source_file))
    write.csv(spectra_input, file = paste(result_dir_spectra, "/spectral_results_for_", file_id, ".csv", sep = ""))
}

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
    no_of_ms2_peaks <- c()
    int <- c() # stores intensity of the MS1 feature
    nx <- 0 # stores number for the ID
    indeX <- 0 # stores number to name the peaklist files

    # pre_mz is a list of precursor m/z
    for (i in pre_mz){


        #filter based on pre mz; sps_all is preprocessed spectra
        sps <- filterPrecursorMzRange(sps_all, i)
        #sps <- filterIntensity(sps, intensity = low_int)

        if (length(sps)>0){

            #ids
            nx <- nx+1
            id_Xx <- paste(file_id,  "M",  as.character(round(i, digits = 0)),
                              "R", as.character(round(median(sps$rtime, na.rm = TRUE), digits = 0)),
                              "ID", as.character(nx), sep = '')
            id_X <- c(id_X, id_Xx)


            #mz
            premz <- c(premz, i)

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


            #peak lists
            # variable for name
            names <- c()

            # create a new directory to store all the peak list txt files
            dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/peakfiles_ms2", sep =""), "."), sep = "")
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
                    fileN1 <- str_replace(fileN, input_dir, ".")
                    ms2Peaks <- c(ms2Peaks, fileN1)
                }
            }
        }
    }



    first_list <- data.frame(cbind(id_X, premz, rtmed, rtmean, int ,col_eng, pol, ms2Peaks))
    write.csv(first_list, file = paste(input_dir, str_remove(paste(result_dir,'/insilico/MS2DATA.csv', sep = ""), "."), sep =""))
    return(first_list)
}



#cam_funcMode <- function(path, pattern = ".mzML"){
    #library("CAMERA")
    # List all files present in QC folder
    #files_QC_N <- list.files(path, pattern = pattern ,full.names=TRUE)

    #for (i in 1:length(files_QC_N)){

        # read each file using Spectra
        #sps_all <- Spectra(files_QC_N[i], backend = MsBackendMzR())


        #if (length(unique(sps_all$polarity)) == 1){
            #if(unique(sps_all$polarity) == 1){
                #Read the same file with MS1 information; note CAMERA reads xcmsSet object
                #xs <- xcmsSet(file = as.character(files_QC_N[i]),
                              #profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                              #mslevel= 1, progressCallback=NULL, polarity="positive",
                              #scanrange = NULL, BPPARAM = bpparam(),
                              #stopOnError = TRUE)
                # Create an xsAnnotate object
                #an <- xsAnnotate(xs)
                # Group based on RT
                #anF <- groupFWHM(an, perfwhm = 0.6)
                # Annotate isotopes
                #anI <- findIsotopes(anF, mzabs = 0.01)
                # Verify grouping
                #anIC <- groupCorr(anI, cor_eic_th = 0.75)
                #Annotate adducts
                #anFA <- findAdducts(anIC, polarity="positive")
                #get a feature list
                #peaklist <- getPeaklist(anFA)
                # add file_origin information
                #peaklist$file_origin <- as.character(files_QC_N[i])
                #file_name <- paste(path,"/posCAMERA_Results_", i, ".csv", sep = "")
                # write individual QC files which are in pos mode (later code will combine them)
                #write.csv(peaklist, file = file_name)
            #}else if (unique(sps_all$polarity) == 0){
                #Read the same file with MS1 information; note CAMERA reads xcmsSet object
                #xs <- xcmsSet(file = as.character(files_QC_N[i]),
                              #profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                              #mslevel= 1, progressCallback=NULL, polarity="negative",
                              #scanrange = NULL, BPPARAM = bpparam(),
                              #stopOnError = TRUE)
                # Create an xsAnnotate object
                #an <- xsAnnotate(xs)
                # Group based on RT
                #anF <- groupFWHM(an, perfwhm = 0.6)
                # Annotate isotopes
                #anI <- findIsotopes(anF, mzabs = 0.01)
                # Verify grouping
                #anIC <- groupCorr(anI, cor_eic_th = 0.75)
                #Annotate adducts
                #anFA <- findAdducts(anIC, polarity="negative")
                #get a feature list
                #peaklist <- getPeaklist(anFA)
                # add file_origin information
                #peaklist$file_origin <- as.character(files_QC_N[i])
                #file_name <- paste(path, "/negCAMERA_Results_", i, ".csv", sep = "")
                # write individual QC files which are in pos mode (later code will combine them)
                #write.csv(peaklist, file = file_name)
            #}
        #}
    #else{
        #pos <- sps_all[sps_all$polarity == 1]
        #neg <- sps_all[sps_all$polarity == 0]
        #file_p <- paste(path, "/QC_280k_pos", i, ".mzML", sep = "")
        #file_n <- paste(path, "/QC_280k_neg", i, ".mzML", sep = "")

        # create new mzML QC files for pos and neg modes from each common QC file
        #export(pos, backend = MsBackendMzR(), file = file_p)
        #export(neg, backend = MsBackendMzR(), file = file_n)
        #Read the same file with MS1 information; note CAMERA reads xcmsSet object
        #xs <- xcmsSet(file = as.character(file_p),
                        #profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                        #mslevel= 1, progressCallback=NULL, polarity="positive",
                        #scanrange = NULL, BPPARAM = bpparam(),
                        #stopOnError = TRUE)
        # Create an xsAnnotate object
        #an <- xsAnnotate(xs)
        # Group based on RT
        #anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes
        #anI <- findIsotopes(anF, mzabs = 0.01)
        # Verify grouping
        #anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts
        #anFA <- findAdducts(anIC, polarity="positive")
        #get a feature list
        #peaklist <- getPeaklist(anFA)
        # add file_origin information
        #peaklist$file_origin <- as.character(files_QC_N[i])
        #file_name <- paste(path, "/posCAMERA_Results_", i, ".csv", sep = "")
        # write individual QC files which are in pos mode (later code will combine them)
        #write.csv(peaklist, file = file_name)

        #Read the same file with MS1 information; note CAMERA reads xcmsSet object
        #xs <- xcmsSet(file = as.character(file_n),
                        #profmethod = "bin", profparam = list(), lockMassFreq=FALSE,
                        #mslevel= 1, progressCallback=NULL, polarity="negative",
                        #scanrange = NULL, BPPARAM = bpparam(),
                        #stopOnError = TRUE)
        # Create an xsAnnotate object
        #an <- xsAnnotate(xs)
        # Group based on RT
        #anF <- groupFWHM(an, perfwhm = 0.6)
        # Annotate isotopes
        #anI <- findIsotopes(anF, mzabs = 0.01)
        # Verify grouping
        #anIC <- groupCorr(anI, cor_eic_th = 0.75)
        #Annotate adducts
        #anFA <- findAdducts(anIC, polarity="negative")
        #get a feature list
        #peaklist <- getPeaklist(anFA)
        # add file_origin information
        #peaklist$file_origin <- as.character(files_QC_N[i])
        #file_name <- paste(path, "/negCAMERA_Results_", i, ".csv", sep = "")
        # write individual QC files which are in pos mode (later code will combine them)
        #write.csv(peaklist, file = file_name)
    #}

    #}

    #detach("package:CAMERA", unload=TRUE)
#}

# ---------- merge_qc ----------

#merge_qc<- function(path){
    # combine all QC which are in positive mode
    #df_pos <- list.files(path, pattern = "posCAMERA_Results_", full.names = TRUE) %>%
        #lapply(read_csv) %>%
        #bind_rows
    # remove any duplicated rows
    #df_pos <- as.data.frame(df_pos[!duplicated(df_pos), ])

    #extract isotope column numbers, the numbers represent the group of isotope
    #nm_p <- regmatches(df_pos[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",df_pos[, "isotopes"]))

    # for all the numbers, extract only first number, since it is the group number,
    # second number can be charge
    #for (i in 1:length(nm_p)){
       # y <- as.numeric(unlist(nm_p[i]))
        #df_pos[i,'istops'] = y[1]
    #}

    # write csv for the combined_camera_pos results
    #write.csv(df_pos, paste(path, "/Combined_Camera_pos.csv", sep = ""))

    # combine all QC which are in negative mode
    #df_neg <- list.files(path, pattern = "negCAMERA_Results_", full.names = TRUE) %>%
        #lapply(read_csv) %>%
        #bind_rows
    # remove any duplicated rows based on mz
    #df_neg <- as.data.frame(df_neg[!duplicated(df_neg), ])

    #extract isotope column numbers, the numbers represent the group of isotope
    #nm_n <- regmatches(df_neg[, "isotopes"],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",df_neg[, "isotopes"]))

    # for all the numbers, extract only first number, since it is the group number,
    # second number can be charge
    #for (i in 1:length(nm_n)){
        #y <- as.numeric(unlist(nm_n[i]))
        #df_neg[i,'istops'] = y[1]
    #}
    # write csv for the combined_camera_neg results
    #write.csv(df_neg, paste(path, "/Combined_Camera_neg.csv", sep = ""))
#}


cam_func <- function(input_dir, f, ms2features){
    modes_file <- read_csv(ms2features)
    mode = unique(modes_file["pol"])
    fl <- paste(input_dir, str_remove(f, "."), sep ="")
    if(mode == "pos"){
        library("CAMERA")
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
        write.csv(peaklist, file = paste(input_dir, str_remove(name, "."), "/CAMERAResults",".csv", sep = ""))
        unloadNamespace("CAMERA")
        unloadNamespace("xcms")
        unloadNamespace("MSnBase")


    }else{
        library("CAMERA")
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
        write.csv(peaklist, file = paste(input_dir, str_remove(name, "."), "/CAMERAResults", ".csv", sep = ""))
        unloadNamespace("CAMERA")
        unloadNamespace("xcms")
        unloadNamespace("MSnBase")
    }
    return(peaklist)
}

# Extract isotopic peaks for each pre_mz
# The input is x = first_list (from ms2peaks function) and y = camera results

ms1_peaks <- function(x, y, result_dir, input_dir, QCfile){
    # store the ms1_peak list path here
    ms1Peaks <- c()
    x = read.csv(x)

    if (QCfile){

        dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/peakfiles_ms1", sep =""), "."), sep = "")
        # create a new directory to store all the peak list txt files
        if (!file.exists(dir_name)){
            dir.create(dir_name, recursive = TRUE)
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

            #if there was only one index
            if (nrow(df_y) == 1){

                # if there was no isotope annotation for that one index
                if (is.na(df_y[1, "istops"])){

                    mz <- df_y[1, "mz"] # save mz
                    int <- df_y[1, "into"] # save intensity
                    no_isotop <- cbind(mz, int) # save as table
                    name_file <- paste(dir_name, "/ms1_peaks_", x[i, 'premz'], "_no_isotopes.txt", sep = "") # save name of the peaklist
                    write.table(no_isotop, name_file, row.names = FALSE, col.names = FALSE) # save peak list
                    name_file1 <- str_replace(name_file, input_dir, ".")
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
                    name_file1 <- str_replace(name_file, input_dir, ".")
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
                    name_file1 <- str_replace(name_file, input_dir, ".")
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
                    name_file1 <- str_replace(name_file, input_dir, ".")
                    ms1Peaks <- c(ms1Peaks, name_file1) # add the path of the peak list to a list
                }
            }
            else if (nrow(df_y)==0){
                ms1Peaks <- c(ms1Peaks, 'no ms1 peaks in QC')
            }
        }
        second_list <- data.frame(cbind(x, ms1Peaks))
        write.csv(second_list, file = paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA.csv', sep = ""), "."), sep =""))
        return(second_list)
    }
    else{

        ms1Peaks <- c(ms1Peaks, 'no ms1 peaks in QC')
        second_list <- data.frame(cbind(x, ms1Peaks))
        write.csv(second_list, file = paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA.csv', sep = ""), "."), sep =""))
        return(second_list)
    }

}


# sirius_postprocess <- function(x, SL = TRUE){
#     feat_scale <- function(p) {
#     (p - min(p)) / (max(p) - min(p))
#     }
#     #the result directory name for each file
#     dir_name <- paste(x, "/insilico/SIRIUS", sep = '')
#     ##result json files in each result directory
#     #for suspect list
#     SL_Param <- list.files(dir_name, pattern = 'SList.json', full.names = TRUE)
#     #for all DB
#     Param <- list.files(dir_name, pattern = 'param.json', full.names = TRUE)
#     #DATA FRAME of both directories for each feature
#     parameter_json <- cbind(SL_Param, Param)
#     #read the msdata csv file that contains all features and their data
#     msdata <- as.data.frame(read.csv(paste(x, "/insilico/MS1DATA.csv", sep = '')))
#     #add new empty columns to store information from SIRIUS results
#     msdata$Adducts <- NA           #adducts
#     msdata$name <- NA              #name of compound
#     msdata$PubChemIDs <- NA        #all pubchem ids
#     msdata$SMILES <- NA            #smiles of compound
#     msdata$Formula <- NA           #formula of compound
#     msdata$FormulaRank <- NA       #formula rank
#     msdata$SIRIUSscore <- NA       #Formula score
#     msdata$CSIFingerIDscore <- NA  #Structure score
#     msdata$SMILESforMCSS <- NA   #SMILES of top scoring candidates; to calculate their tanimoto later
#     msdata$exp_int <- NA           #explained intensity of formula
#     msdata$dir <- NA               #directory for results either the strcuture or formula candidates
#     msdata$Result <- NA             #the type of candidate, either from Suspect list, strcuture command or formula command
#     # for all entries in the json result files for each each feature
#     for (i in 1:nrow(parameter_json)){

#         # find the number of the row corresponding to the row in msdata that has the same precursor m/z as json result files
#         rowMS <- msdata[grepl(str_match(as.character(parameter_json[i, 'Param']), "MS1p_\\s*(.*?)\\s*_SIRIUS")[2] ,msdata$premz), ]

#         if (SL){

#             # file path for strcuture candidate from Suspect List json folder
#             str_canS <- paste(list.dirs(parameter_json[i,'SL_Param'])[2], '/structure_candidates.tsv', sep = '')
#             # file path for formula candidate from Param json folder
#             for_canS <- paste(list.dirs(parameter_json[i,'SL_Param'])[2], '/formula_candidates.tsv', sep = '')
#             # file path for strcuture candidate from Suspect List json folder
#             str_can <- paste(list.dirs(parameter_json[i,'Param'])[2], '/structure_candidates.tsv', sep = '')
#             # file path for formula candidate from Param json folder
#             for_can <- paste(list.dirs(parameter_json[i,'Param'])[2], '/formula_candidates.tsv', sep = '')

#             # if the strcuture candidate file exists
#             if (file.exists(str_canS) && file.exists(str_can)){

#                 # read the corresponding structure and formula candidate files
#                 str_canSL <- as.data.frame(read_tsv(str_canS))
#                 for_canSL <- as.data.frame(read_tsv(for_canS))

#                 # read the corresponding structure and formula candidate files
#                 str_canP <- as.data.frame(read_tsv(str_can))
#                 for_canP <- as.data.frame(read_tsv(for_can))

#                 # if the strcuture candidate file contains 1 or more rows, it has detected a candidate from suspect list, add relevant info

#                 if (nrow(str_canSL) >= 1){

#                     if (str_canSL[1, 'CSI:FingerIDScore'] > str_canP[1, 'CSI:FingerIDScore']){
#                         # information from structure candidate file
#                         msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canSL[1, 'adduct']
#                         msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canSL[1, 'name']
#                         msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canSL[1, 'pubchemids']
#                         msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canSL[1, 'smiles']
#                         msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canSL[1, 'molecularFormula']
#                         msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canSL[1, 'formulaRank']
#                         msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canSL[1, 'CSI:FingerIDScore']

#                         # information from formula candidate file
#                         formulaRow <- which(for_canSL[,'rank'] == str_canSL[1, 'formulaRank'])
#                         msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canSL[formulaRow, 'SiriusScore']
#                         msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canSL[formulaRow, 'explainedIntensity']

#                         # other info
#                         msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_SL'
#                         msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_canS
#                     }
#                     else{

#                         # if the structure candidate file contains 1 row, it has detected a candidate from all DBs, add relevant info
#                         if (nrow(str_canP) == 1){

#                             # information from structure candidate file
#                             msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP[1, 'adduct']
#                             msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP[1, 'name']
#                             msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP[1, 'pubchemids']
#                             msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP[1, 'smiles']
#                             msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP[1, 'molecularFormula']
#                             msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP[1, 'formulaRank']
#                             msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP[1, 'CSI:FingerIDScore']

#                             # information from formula candidate file
#                             formulaRow1 <- which(for_canP[,'rank'] == str_canP[1, 'formulaRank'])
#                             msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[formulaRow1, 'SiriusScore']
#                             msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[formulaRow1, 'explainedIntensity']

#                             # other info
#                             msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                             msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can
#                         }
#                         # if the structure candidate file contains more rows, extract SMILES of top candidates and check their similarity later
#                         else if (nrow(str_canP) > 1){

#                             # information from structure candidate file
#                             msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP[1, 'adduct']
#                             msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP[1, 'name']
#                             msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP[1, 'pubchemids']
#                             msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP[1, 'smiles']
#                             msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP[1, 'molecularFormula']
#                             msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP[1, 'formulaRank']
#                             msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP[1, 'CSI:FingerIDScore']

#                             # information from formula candidate file, take info from the formula rank that corresponds to the formula rank with the top strcuture candidate
#                             formulaRow2 <- which(for_canP[,'rank'] == str_canP[1, 'formulaRank'])
#                             msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[formulaRow2, 'SiriusScore']
#                             msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[formulaRow2, 'explainedIntensity']

#                             # other info
#                             msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                             msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can

#                             # normalize the CSI:FingerIDScores
#                             norm_score <- feat_scale(str_canP[,"CSI:FingerIDScore"])
#                             # store the upper quartile
#                             upper_quartile <- str_canP[which(norm_score > as.numeric(quantile(norm_score)[4])), "smiles"]

#                             # if the upper quartile has more than 5 candidates, then just take the top 5 candidates
#                             if (length(upper_quartile) > 5){
#                                 upper_quartile <- upper_quartile[1:5]
#                             }
#                             # save the top candidates SMILES, to check similarity later with rdkit in Python
#                             msdata[as.numeric(rownames(rowMS)), 'SMILESforMCSS'] <- paste(upper_quartile, collapse = '|')

#                         }
#                         # if the structure candidate file is empty, take information from just the formula candidate file
#                         else if (nrow(str_canP) == 0){

#                             # if formula candidate file is not empty
#                             if (nrow(for_canP) >= 1){

#                                 # information from formula candidate file
#                                 msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP[1, 'adduct']
#                                 msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP[1, 'molecularFormula']
#                                 msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP[1, 'rank']
#                                 msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[1, 'explainedIntensity']
#                                 msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[1, 'SiriusScore']
#                                 # other info
#                                 msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                                 msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can

#                             }
#                         }

#                         # if the structure candidate from all DBs does not exist
#                         else{
#                             # check if the formula candidate file exists
#                             if (file.exists(for_can)){
#                                 for_canF1 <- as.data.frame(read_tsv(for_can))

#                                 # if formula candidate file is not empty
#                                 if (nrow(for_canF1)>= 1){

#                                     # information from formula candidate file
#                                     msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canF1[1, 'adduct']
#                                     msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canF1[1, 'molecularFormula']
#                                     msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canF1[1, 'rank']
#                                     msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canF1[1, 'explainedIntensity']
#                                     msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canF1[1, 'SiriusScore']
#                                     # other info
#                                     msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                                     msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can
#                                 }
#                             }
#                         }
#                     }
#                 }
#                 # if it's empty, move onto the All DB result folder called PARAM here
#                 else{

#                     # if the structure candidate file contains 1 row, it has detected a candidate from all DBs, add relevant info
#                     if (nrow(str_canP) == 1){

#                         # information from structure candidate file
#                         msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP[1, 'adduct']
#                         msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP[1, 'name']
#                         msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP[1, 'pubchemids']
#                         msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP[1, 'smiles']
#                         msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP[1, 'molecularFormula']
#                         msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP[1, 'formulaRank']
#                         msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP[1, 'CSI:FingerIDScore']

#                         # information from formula candidate file
#                         formulaRow1 <- which(for_canP[,'rank'] == str_canP[1, 'formulaRank'])
#                         msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[formulaRow1, 'SiriusScore']
#                         msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[formulaRow1, 'explainedIntensity']

#                         # other info
#                         msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                         msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can
#                     }
#                     # if the structure candidate file contains more rows, extract SMILES of top candidates and check their similarity later
#                     else if (nrow(str_canP) > 1){

#                         # information from structure candidate file
#                         msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP[1, 'adduct']
#                         msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP[1, 'name']
#                         msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP[1, 'pubchemids']
#                         msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP[1, 'smiles']
#                         msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP[1, 'molecularFormula']
#                         msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP[1, 'formulaRank']
#                         msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP[1, 'CSI:FingerIDScore']

#                         # information from formula candidate file, take info from the formula rank that corresponds to the formula rank with the top strcuture candidate
#                         formulaRow2 <- which(for_canP[,'rank'] == str_canP[1, 'formulaRank'])
#                         msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[formulaRow2, 'SiriusScore']
#                         msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[formulaRow2, 'explainedIntensity']

#                         # other info
#                         msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                         msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can

#                         # normalize the CSI:FingerIDScores
#                         norm_score <- feat_scale(str_canP[,"CSI:FingerIDScore"])
#                         # store the upper quartile
#                         upper_quartile <- str_canP[which(norm_score > as.numeric(quantile(norm_score)[4])), "smiles"]

#                         # if the upper quartile has more than 5 candidates, then just take the top 5 candidates
#                         if (length(upper_quartile) > 5){
#                             upper_quartile <- upper_quartile[1:5]
#                         }
#                         # save the top candidates SMILES, to check similarity later with rdkit in Python
#                         msdata[as.numeric(rownames(rowMS)), 'SMILESforMCSS'] <- paste(upper_quartile, collapse = '|')

#                     }
#                     # if the structure candidate file is empty, take information from just the formula candidate file
#                     else if (nrow(str_canP) == 0){

#                         # if formula candidate file is not empty
#                         if (nrow(for_canP) >= 1){

#                             # information from formula candidate file
#                             msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP[1, 'adduct']
#                             msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP[1, 'molecularFormula']
#                             msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP[1, 'rank']
#                             msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP[1, 'explainedIntensity']
#                             msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP[1, 'SiriusScore']
#                             # other info
#                             msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                             msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can

#                         }
#                     }

#                     # if the structure candidate from all DBs does not exist
#                     else{
#                         # check if the formula candidate file exists
#                         if (file.exists(for_can)){
#                             for_canF1 <- as.data.frame(read_tsv(for_can))

#                             # if formula candidate file is not empty
#                             if (nrow(for_canF1)>= 1){

#                                 # information from formula candidate file
#                                 msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canF1[1, 'adduct']
#                                 msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canF1[1, 'molecularFormula']
#                                 msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canF1[1, 'rank']
#                                 msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canF1[1, 'explainedIntensity']
#                                 msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canF1[1, 'SiriusScore']
#                                 # other info
#                                 msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                                 msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can
#                             }
#                         }
#                     }
#                 }
#             }
#             else{
#                 # directory for structure candidate file
#                 str_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/structure_candidates.tsv', sep = '')
#                 # directory for formula candidate file
#                 for_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/formula_candidates.tsv', sep = '')

#                 # if the structure candidate file exists
#                 if (file.exists(str_can1)){

#                     # read the structure file from All Dbs (PARAM) json file
#                     str_canP1 <- as.data.frame(read_tsv(str_can1))
#                     # read the formula file from All Dbs (PARAM) json file
#                     for_canP1 <- as.data.frame(read_tsv(for_can1))

#                     #if the structure candidate file has one candidate, it has detected a candidate from all DBs, add relevant info
#                     if (nrow(str_canP1) == 1){

#                         # information from structure candidate file
#                         msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
#                         msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
#                         msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
#                         msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
#                         msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
#                         msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
#                         msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
#                         msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
#                         # information from formula candidate file
#                         formulaRow3 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
#                         msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow3, 'SiriusScore']
#                         msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow3, 'explainedIntensity']
#                         # other info
#                         msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                         msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1

#                     }
#                     # if the strcuture cabdidate file has more than 1 candidates, it has detected candidates from all DBs, add relevant info
#                     else if (nrow(str_canP1) > 1){

#                         # information from structure candidate file
#                         msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
#                         msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
#                         msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
#                         msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
#                         msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
#                         msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
#                         msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
#                         # information from formula candidate file
#                         formulaRow4 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
#                         msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow4, 'SiriusScore']
#                         msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow4, 'explainedIntensity']
#                         # other info
#                         msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                         msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1

#                         # normalize the CSI:FingerIDScores
#                         norm_score1 <- feat_scale(str_canP1[,"CSI:FingerIDScore"])
#                         # store the upper quartile
#                         upper_quartile1 <- str_canP1[which(norm_score1 > as.numeric(quantile(norm_score1)[4])), "smiles"]
#                         # if the upper quartile has more than 5 candidates, then just take the top 5 candidates
#                         if (length(upper_quartile1) > 5){
#                             upper_quartile1 <- upper_quartile1[1:5]
#                         }
#                         # save the top candidates SMILES, to check similarity later with rdkit in Python
#                         msdata[as.numeric(rownames(rowMS)), 'SMILESforMCSS'] <- paste(upper_quartile1, collapse = '|')
#                     }
#                     #
#                     else if (nrow(str_canP1) == 0){
#                         if (file.exists(for_can1)){
#                             for_canP2 <- as.data.frame(read_tsv(for_can1))
#                             if (nrow(for_canP2)>= 1){

#                                 # information from formula candidate file
#                                 msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP2[1, 'adduct']
#                                 msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP2[1, 'molecularFormula']
#                                 msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP2[1, 'rank']
#                                 msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP2[1, 'explainedIntensity']
#                                 msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP2[1, 'SiriusScore']
#                                 # other info
#                                 msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                                 msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
#                             }
#                         }
#                     }
#                 }

#                 # if the structure candidate file doesn't exists (and no str_sl exists)
#                 else{
#                     # if formula candidate file exists
#                     if (file.exists(for_can1)){
#                         for_canP3 <- as.data.frame(read_tsv(for_can1))
#                         if (nrow(for_canP3) >= 1){
#                             # information from formula candidate file
#                             msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP3[1, 'adduct']
#                             msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP3[1, 'molecularFormula']
#                             msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP3[1, 'rank']
#                             msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP3[1, 'explainedIntensity']
#                             msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP3[1, 'SiriusScore']
#                             # other info
#                             msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                             msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
#                         }
#                     }
#                 }
#             }
#         }# SL IS false
#         else{
#             # directory for structure candidate file
#             str_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/structure_candidates.tsv', sep = '')
#             # directory for formula candidate file
#             for_can1 <- paste(list.dirs(parameter_json[i,'Param'])[2], '/formula_candidates.tsv', sep = '')

#             # if the structure candidate file exists
#             if (file.exists(str_can1)){

#                 # read the structure file from All Dbs (PARAM) json file
#                 str_canP1 <- as.data.frame(read_tsv(str_can1))
#                 # read the formula file from All Dbs (PARAM) json file
#                 for_canP1 <- as.data.frame(read_tsv(for_can1))

#                 #if the structure candidate file has one candidate, it has detected a candidate from all DBs, add relevant info
#                 if (nrow(str_canP1) == 1){

#                     # information from structure candidate file
#                     msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
#                     msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
#                     msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
#                     msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
#                     msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
#                     msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
#                     msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
#                     msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
#                     # information from formula candidate file
#                     formulaRow3 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
#                     msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow3, 'SiriusScore']
#                     msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow3, 'explainedIntensity']
#                     # other info
#                     msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                     msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1

#                 }
#                 # if the strcuture cabdidate file has more than 1 candidates, it has detected candidates from all DBs, add relevant info
#                 else if (nrow(str_canP1) > 1){

#                     # information from structure candidate file
#                     msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- str_canP1[1, 'adduct']
#                     msdata[as.numeric(rownames(rowMS)), 'name'] <- str_canP1[1, 'name']
#                     msdata[as.numeric(rownames(rowMS)), 'PubChemIDs'] <- str_canP1[1, 'pubchemids']
#                     msdata[as.numeric(rownames(rowMS)), 'SMILES'] <- str_canP1[1, 'smiles']
#                     msdata[as.numeric(rownames(rowMS)), 'Formula'] <- str_canP1[1, 'molecularFormula']
#                     msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- str_canP1[1, 'formulaRank']
#                     msdata[as.numeric(rownames(rowMS)), 'CSIFingerIDscore'] <- str_canP1[1, 'CSI:FingerIDScore']
#                     # information from formula candidate file
#                     formulaRow4 <- which(for_canP1[,'rank'] == str_canP1[1, 'formulaRank'])
#                     msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP1[formulaRow4, 'SiriusScore']
#                     msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP1[formulaRow4, 'explainedIntensity']
#                     # other info
#                     msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_STR'
#                     msdata[as.numeric(rownames(rowMS)), 'dir'] <- str_can1

#                     # normalize the CSI:FingerIDScores
#                     norm_score1 <- feat_scale(str_canP1[,"CSI:FingerIDScore"])
#                     # store the upper quartile
#                     upper_quartile1 <- str_canP1[which(norm_score1 > as.numeric(quantile(norm_score1)[4])), "smiles"]
#                     # if the upper quartile has more than 5 candidates, then just take the top 5 candidates
#                     if (length(upper_quartile1) > 5){
#                         upper_quartile1 <- upper_quartile1[1:5]
#                     }
#                     # save the top candidates SMILES, to check similarity later with rdkit in Python
#                     msdata[as.numeric(rownames(rowMS)), 'SMILESforMCSS'] <- paste(upper_quartile1, collapse = '|')
#                 }
#                 #
#                 else if (nrow(str_canP1) == 0){
#                     if (file.exists(for_can1)){
#                         for_canP2 <- as.data.frame(read_tsv(for_can1))
#                         if (nrow(for_canP2)>= 1){

#                             # information from formula candidate file
#                             msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP2[1, 'adduct']
#                             msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP2[1, 'molecularFormula']
#                             msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP2[1, 'rank']
#                             msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP2[1, 'explainedIntensity']
#                             msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP2[1, 'SiriusScore']
#                             # other info
#                             msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                             msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
#                         }
#                     }
#                 }
#             }

#             # if the structure candidate file doesn't exists (and no str_sl exists)
#             else{
#                 # if formula candidate file exists
#                 if (file.exists(for_can1)){
#                     for_canP3 <- as.data.frame(read_tsv(for_can1))
#                     if (nrow(for_canP3) >= 1){
#                         # information from formula candidate file
#                         msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP3[1, 'adduct']
#                         msdata[as.numeric(rownames(rowMS)), 'Formula'] <- for_canP3[1, 'molecularFormula']
#                         msdata[as.numeric(rownames(rowMS)), 'FormulaRank'] <- for_canP3[1, 'rank']
#                         msdata[as.numeric(rownames(rowMS)), 'exp_int'] <- for_canP3[1, 'explainedIntensity']
#                         msdata[as.numeric(rownames(rowMS)), 'SIRIUSscore'] <- for_canP3[1, 'SiriusScore']
#                         # other info
#                         msdata[as.numeric(rownames(rowMS)), 'Result'] <- 'SIRIUS_FOR'
#                         msdata[as.numeric(rownames(rowMS)), 'dir'] <- for_can1
#                     }
#                 }
#             }
#         }
#     }
#     write.csv(msdata, paste(x, "/insilico/MS1DATAsirius.csv", sep = ''))
#     return(msdata)
# }




# input x is result from either ms1_peaks
# SL is if a suspect list is present

sirius_param <- function(x, result_dir, input_dir, SL = FALSE, collision_info = FALSE){

    dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/SIRIUS", sep =""), "."), sep = "")
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
            fileSR <- paste(str_sub(fileR, end=-4),'.json', sep = '')
            fileSRS <- paste(str_sub(fileR, end=-4), 'SList.json', sep = '')
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
#             writeLines(">ms1",con=file.conn)
#             writeLines(paste(x[i,"premz"], x[i,"int"] ,sep=" "),con=file.conn)

            #ms2
            if (collision_info){
                writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            }else{
                writeLines(">ms2" ,con=file.conn)
            }
            
            ms2pk <- paste(input_dir, str_remove(x[i,"ms2Peaks"], "."), sep ="")

            peak<- read.table(ms2pk)
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
            fileSR <- paste(str_sub(fileR, end=-4),'.json', sep = '')
            fileSRS <- paste(str_sub(fileR, end=-4), 'SList.json', sep = '')
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

            ms1pk <- paste(input_dir, str_remove(x[i,"ms1Peaks"], "."), sep ="")
            peakms1<- read.table(ms1pk)

            for (l in 1:length(peakms1[,1])){
                writeLines(paste(as.character(peakms1[l,1]),as.character(peakms1[l,2]), sep =" "), con=file.conn)
            }

            #ms2
            if (collision_info){
                writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            }else{
                writeLines(">ms2" ,con=file.conn)
            }

            ms2pk <- paste(input_dir, str_remove(x[i,"ms2Peaks"], "."), sep ="")

            peakms2<- read.table(ms2pk)

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
            fileSR <- paste(str_sub(fileR, end=-4),'.json', sep = '')
            fileSRS <- paste(str_sub(fileR, end=-4), 'SList.json', sep = '')
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

            ms1pk <- paste(input_dir, str_remove(x[i,"ms1Peaks"], "."), sep ="")
            peakms1<- read.table(ms1pk)

            for (l in 1:length(peakms1[,1])){
                writeLines(paste(as.character(peakms1[l,1]),as.character(peakms1[l,2]), sep =" "), con=file.conn)
            }

            #ms2
            if (collision_info){
                writeLines(paste(">collision", paste(x[i,"col_eng"],"eV", sep =''),sep=" "),con=file.conn)
            }else{
                writeLines(">ms2" ,con=file.conn)
            }

            ms2pk <- paste(input_dir, str_remove(x[i,"ms2Peaks"], "."), sep ="")

            peakms2<- read.table(ms2pk)

            for (k in 1:length(peakms2[,1])){
                writeLines(paste(as.character(peakms2[k,1]),as.character(peakms2[k,2]), sep =" "), con=file.conn)
            }

            close(file.conn)
            parameter_file <- c(parameter_file,file.conn)

        }
    }
    if (SL){

        in_out_file <- data.frame(cbind(sirius_param_file, outputNames, outputNamesSL, isotopes))

        write.table(in_out_file, paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA_SiriusPandSL.tsv', sep = ""), "."), sep =""), sep = "\t")
        return(in_out_file)

    }
    else{
        in_out_file <- data.frame(cbind(sirius_param_file, outputNames, isotopes))

        write.table(in_out_file, paste(input_dir, str_remove(paste(result_dir,'/insilico/MS1DATA_SiriusP.tsv', sep = ""), "."), sep =""), sep = "\t")
        return(in_out_file)

    }

}


run_sirius <- function(files, ppm_max = 5, ppm_max_ms2 = 15, QC = TRUE, SL = TRUE, SL_path, candidates = 30, profile, db){
    files <- read.csv(files, sep = "\t")

    for (b in 1:nrow(files)){

        if (QC){
            if (is.na(files[b, "isotopes"])){
                command1 = c("sirius",  "--input", files[b, "sirius_param_file"], 
                "--output", files[b, "outputNames"],
                "formula", "--profile", profile, "--no-isotope-filter",  "--no-isotope-score", "--candidates", 
                candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2, "structure", "--database", 
                db ,"canopus")
                print(command1)
                print("command1")
                system2(command1)
                if(!file.exists(files[b, "outputNames"])){
                    command2 = c("sirius", "--input", files[b, "sirius_param_file"], "--output", 
                    files[b, "outputNames"],
                    "formula", "--profile", 
                    profile, "--no-isotope-filter", "--no-isotope-score", "--candidates", 
                    candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"structure",  "--database", 
                    db ,"canopus")
                    print(command2)
                    print("command2")
                    system2(command2)
                }
                if(SL){
                    command3 = c("sirius", "--input", files[b, "sirius_param_file"], 
                    "--output", files[b, "outputNamesSL"],
                    "formula", "--profile", profile, "--no-isotope-filter", "--no-isotope-score", "--candidates", 
                    candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,"structure", "--database",  
                    SL_path ,"canopus")
                    
                    print(command3)
                    print("command3")
                    system2(command3)
                }

            }
            else if(files[b, "isotopes"] == "present"){

                command4 = c("sirius", "--input", files[b, "sirius_param_file"], "--output", 
                files[b, "outputNames"],
                "formula", "--profile", profile, "--candidates", candidates, 
                "--ppm-max", ppm_max,"--ppm-max-ms2", ppm_max_ms2,"structure", "--database", db ,"canopus"
                            )
                print(command4)
                print("command4")
                system2(command4)
                if(!file.exists(files[b, "outputNames"])){
                    command5 = c("sirius",  "--input", files[b, "sirius_param_file"], 
                    "--output", files[b, "outputNames"],
                    "formula", "--profile", profile, "--candidates", candidates, 
                    "--ppm-max", ppm_max,"--ppm-max-ms2", ppm_max_ms2, "structure", "--database", 
                    db ,"canopus")
                    print(command5)
                    print("command5")
                    system2(command5)
                }
                if(SL){
                    command6 = c("sirius", "--input", files[b, "sirius_param_file"], 
                    "--output", 
                    files[b, "outputNamesSL"],
                    "formula", "--profile", 
                    profile, "--candidates", 
                    candidates, "--ppm-max", 
                    ppm_max, "--ppm-max-ms2",  
                    ppm_max_ms2, "structure",  "--database",  SL_path ,"canopus")
                    print(command6)
                    print("command6")
                    system2(command6)
                }
            }
        }
        else{
            command7 = c("sirius", "--input", files[b, "sirius_param_file"], "--output", files[b, "outputNames"],
                        "formula", "--profile", profile, "--no-isotope-filter",
                        "--no-isotope-score",
                        "--candidates", 
                        candidates, 
                        "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,
                        "structure", "--database", db ,"canopus")
            print(command7)
            print("command7")
            system2(command7)
            if(!file.exists(files[b, "outputNames"])){
                command8 = c("sirius", 
                "--input", 
                files[b, "sirius_param_file"],
                "--output", 
                files[b, "outputNames"],
                "formula",
                "--profile", 
                profile, "--no-isotope-filter", "--no-isotope-score", "--candidates", 
                candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,
                "structure", "--database", db ,"canopus")
                print(command8)
                print("command8")
                system2(command8)
            }
            if(SL){
                command9 = c("sirius", "--input", files[b, "sirius_param_file"], 
                "--output", files[b, "outputNamesSL"], "formula", "--profile", 
                profile, "--no-isotope-filter", "--no-isotope-score", "--candidates", 
                candidates, "--ppm-max", ppm_max, "--ppm-max-ms2",  ppm_max_ms2,
                "structure", "--database ",  SL_path ,"canopus")
                
                print(command9)
                print("command9")
                system2(command9)
            }
        }
        Sys.sleep(5)

    }


}




sirius_adduct <- function(input_dir, x, SL = TRUE){

    feat_scale <- function(p) {
    (p - min(p)) / (max(p) - min(p))
    }

    #the result directory name for each file
    dir_name <- paste(input_dir, str_remove(x, "."),'/insilico/SIRIUS', sep ="")

    #for all DB
    Param <- list.files(dir_name, pattern = 'param.json', full.names = TRUE)
    #read the msdata csv file that contains all features and their data
    msdatacsv <- read.csv(paste(input_dir, str_remove(x, "."), "/insilico/MS1DATA.csv", sep = ''))
    msdata <- as.data.frame(msdatacsv)
    #add new empty columns to store information from SIRIUS results
    msdata$Adducts <- NA           #adducts

    # for all entries in the json result files for each each feature
    for (i in Param){

        # find the number of the row corresponding to the row in msdata that has the same precursor m/z as json result files
        rowMS <- msdata[grepl(str_match(as.character(i), "MS1p_\\s*(.*?)\\s*_SIRIUS")[2] ,msdata$premz), ]

        # directory for formula candidate file
        for_can1 <- paste(list.dirs(i)[2], '/formula_candidates.tsv', sep = '')
        if (file.exists(for_can1)){
            # read the formula file from All Dbs (PARAM) json file
            for_canP1 <- as.data.frame(read_tsv(for_can1))

            if (nrow(for_canP1)>= 1){

                # information from formula candidate file
                msdata[as.numeric(rownames(rowMS)), 'Adducts'] <- for_canP1[1, 'adduct']
            }
        }
    }
    write.csv(msdata, paste(input_dir, str_remove(x, "."), "/insilico/MS1DATAsirius.csv", sep = ''))
    return(msdata)
}



metfrag_param <- function(x, result_dir, input_dir, sl_mtfrag, SL = TRUE, ppm_max = 5, ppm_max_ms2= 15){



    if (file.exists(paste(input_dir, "/MetFrag_AdductTypes.csv", sep = ""))){
        adducts <- paste(input_dir, "/MetFrag_AdductTypes.csv", sep = "")
    }
    else{
        system(paste("wget -P",
                    input_dir,
                    "https://github.com/schymane/ReSOLUTION/raw/master/inst/extdata/MetFrag_AdductTypes.csv",
                    sep = " "))
        adducts <- paste(input_dir, "/MetFrag_AdductTypes.csv", sep = "")
    }

    x <- read.csv(x)

    dir_name <- paste(input_dir, str_remove(paste(result_dir, "/insilico/MetFrag", sep =""), "."), sep = "")

    if (!file.exists(dir_name)){
        dir.create(dir_name, recursive = TRUE) ##create folder
    }

    db <- "PubChem"

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


                peakspath <- str_replace(x[j, "ms2Peaks"], ".", input_dir)


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
                resultspath <- str_replace(result_dir, ".", input_dir)
                writeLines(paste("ResultsPath = ", resultspath, "/insilico/MetFrag/", sep = ''),con=file.conn)

                writeLines("MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter",con=file.conn)
                writeLines("MetFragPostProcessingCandidateFilter = InChIKeyFilter",con=file.conn)
                writeLines("#MaximumTreeDepth = 2",con=file.conn)
                writeLines("#NumberThreads = 1",con=file.conn)

                close(file.conn)
                parameter_file <- c(parameter_file,file.conn)

            }
        }
    }

    write.table(metfrag_param_file, file = paste(input_dir, str_remove(paste(result_dir, "/insilico/metparam_list.txt", sep =""), "."), sep = ""), sep = "/t", row.names = FALSE, col.names = FALSE)
    return(metfrag_param_file)
}

run_metfrag <- function(met_param){
    #line <- system.time({

        if (file.exists(paste(input_dir, "/MetFragCommandLine-2.5.0.jar", sep = ""))){
            MetFragjarFile <- paste(input_dir, "/MetFragCommandLine-2.5.0.jar", sep = "")
        }
        else{
            system(paste("wget -P",
                        input_dir,
                        "https://github.com/ipb-halle/MetFragRelaunched/releases/download/v.2.5.0/MetFragCommandLine-2.5.0.jar",
                        sep = " "))
            MetFragjarFile <- paste(input_dir, "/MetFragCommandLine-2.5.0.jar", sep = "")
        }

        start_time <- Sys.time()

        filesmet_param <- read.table(met_param)

        for (files in filesmet_param[[1]]){
            system(paste("java -jar",  MetFragjarFile , files))
            Sys.sleep(5)
        }
    #})

#     if (file.exists(paste(input_dir, "/summaryFile.txt", sep = ""))){

#         line1 = paste("Compound Database Dereplication with MetFrag took", line[1], "of seconds.", sep = " ")
#         write(line1 ,file=paste(input_dir, "/summaryFile.txt", sep = ""),append=TRUE)

#     }
#     else{
#         fileConn<-file(paste(input_dir, "/summaryFile.txt", sep = ""))

#         writeLines(paste("Compound Database Dereplication with Metfrag took", line[1], "of seconds.", sep = " "), fileConn)

#         close(fileConn)
#     }

}



library(parallel)
library(doParallel)
library(future)
library(iterators)
library(listenv)

options(future.globals.maxSize = 8 * 1024^3) # increase dataset size limit taken by future to 8GB

# detects number of cores
n.cores <- parallel::detectCores()


plan(list(
  tweak(multisession, workers = ((n.cores + 5) %/% 3) %/% 2),
  tweak(multisession, workers = ((n.cores + 5) %/% 3) %/% 2),
  tweak(multisession, workers = 3)
))

# Start time
start.time <- Sys.time()

args = commandArgs(trailingOnly = TRUE)



# input directory
input_dir <- args[1]
#input_dir


# generate a table for all the input .mzML files, along with creating a result directory for all each.mzML file
input_table <- data.frame(ms2_rfilename(input_dir))
#input_table

# take the list of .mzML files as an environment and parallelizes the computation
# for each .mzML file in the for loop below
input_table_idxs <- listenv()

# for loop, to compute MAW-R for all given input .mzML files
for (i in 1:nrow(input_table)){

    # for each index, taken as environment using the listenv library
    input_table_idxs[[i]] <- future({
        
        ##Preprocess and Read the mzML files
        spec_pr <- spec_Processing(input_dir,
                                   input_table[i, "mzml_files"],
                                   input_table[i, "ResultFileNames"])

        #perform dereplication with all dbs (gnps, massbank, and hmdb)
        df_derep <- spec_dereplication_file(mzml_file = input_table[i, "mzml_files"],
                                            pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "."), sep =""),
                                            proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "."), sep =""),
                                            db = "all",
                                            result_dir = input_table[i, "ResultFileNames"],
                                            file_id = input_table[i, "File_id"],
                                            input_dir,
                                            no_of_candidates = 50,
                                            ppmx = 15)
        
        #Extract MS2 peak lists
        spec_pr2 <- ms2_peaks(pre_tbl = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/premz_list.txt", sep = ""), "."), sep =""),
                              proc_mzml = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/processedSpectra.mzML", sep = ""), "."), sep =""),
                              input_dir,
                              result_dir = input_table[i, "ResultFileNames"],
                             file_id = input_table[i, "File_id"])
        

        # camera results for isotopes
        # cam_res <- cam_func(input_dir,
                            # f = input_table[i, "mzml_files"], 
                            # ms2features = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"], "/insilico/MS2DATA.csv", sep = ""), "."), sep = ""))
        

        # Extract MS1 peaks or isotopic peaks
        ms1p <- ms1_peaks(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS2DATA.csv', sep = ""), "."), sep =""),
                          y = NA, 
                          input_table[i, "ResultFileNames"],
                          input_dir,
                          QCfile = FALSE)


        #prepare sirius parameter files
        sirius_param_files <- sirius_param(x = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA.csv', sep = ""), "."), sep =""),
                                           result_dir = input_table[i, 'ResultFileNames'],
                                           input_dir,
                                           SL = FALSE)
        # Run sirius
        #run_sirius(files = paste(input_dir, str_remove(paste(input_table[i, "ResultFileNames"],'/insilico/MS1DATA_SiriusP.tsv', sep = ""), "."), sep =""),
        #           ppm_max = 5,
        #           ppm_max_ms2 = 15,
        #           QC = FALSE,
        #           SL = FALSE,
        #           SL_path = NA,
        #           candidates = 30,
        #          profile = "orbitrap", 
        #          db = "ALL")
   }) #end input_table_idxs future
}



# Run the workflow
input_table_idxs <- as.list(input_table_idxs)
v_input_table_idxs <- future::value(input_table_idxs)

# End time
end.time <- Sys.time()

# Time taken to run the analysis for MAW-R
time.taken <- end.time - start.time
print(time.taken)




