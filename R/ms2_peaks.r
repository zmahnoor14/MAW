#! /usr/bin/Rscript

#' @title Extract MS2 Fragment peaks and save as text files
#'
#' @description
#'
#' For each precursor m/z in the input .mzML file, filter all spectra  
#' with that precursor m/z. Extract retention time, polarity, collision 
#' energy, and maximum intensity from these spectra. Create a directory 
#' named result_dir/insilico/peakfiles_ms2/ to store the ms2 peak list. 
#' To generate these peaklists, unlist the MS2 peaks from each spectrum 
#' and combine these peaks using the function combinePeaks from Spectra 
#' package. Give each peak list a number and add it to the file name, 
#' where it is stored, such as peaks_01.txt. At the end store all this 
#' information into a dataframe and return the dataframe
#'

#'
#' @param pre_tbl txt file of list of precursor m/z (from func spec_Processing)
#'
#' @param proc_mzml processed mzML spectra file (from func spec_Processing)
#'
#' @param input_dir is input directory with all input files
#' 
#' @param result_dir result_dir is the result directory for an input file


#' @return
#' 
#' Peaklists for each precursor m/z, stored in a txt file.
#' A dataframe consisting of features from the input files

#'
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' ms2_peaks("/usr/project/file1/premz_list.txt", "/usr/project/file1/processedSpectra.mzmL",
#' "usr/project/", 'usr/project/file_09/')

# ---------- Preparations ----------
# Load libraries
library("Spectra")
library("stringr")

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)

pre_tbl <- as.character(args[1])
proc_mzml <- as.character(args[2])
input_dir <- as.character(args[3])
result_dir <- as.character(args[4])
file_id <- as.character(args[5])

# ---------- ms2_peaks ----------


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
        id_Xx <- paste(file_id,  "M",  as.character(round(i, digits = 0)), 
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

ms2_peaks(pre_tbl, proc_mzml, input_dir, result_dir)

