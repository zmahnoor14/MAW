###---- library ----
library(RColorBrewer)           # For colors
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(vegan)                  # For shannon diversity
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(circlize)               # For sunburst plot
library(plotrix)                # For sunburst plot
library(caret)                  # Swiss-army knife for statistics
library(pROC)                   # Evaluation metrics
library(PRROC)                  # Evaluation metrics
library(multiROC)               # Evaluation metrics
library(plotly)                 # For creating html plots
library(htmlwidgets)            # For creating html plots
library(stringr)              
source("https://raw.githubusercontent.com/ipb-halle/iESTIMATE/main/R/_functions.r")


# ---- data preparation ----
data_preparation <- function(files, phenodata, result_dir_name, plots_dir_name, chrom_run_sec, msLevel = NULL) {
  # files = list of mzML files
  # phenodata = csv file of phenodata, including (1) sample_name, (2) sample_group, (3) sample_description
  # result_dir_name = name for directory to store all result files and objects
  # plots_dir_name = name for directory to store all plots 
  # chrom_run_sec = length of the chromatography run in seconds
  # msLevel = level of MS for files, default is ms level 1
  
  # create results directory
  if (dir.exists(paste(getwd(), result_dir_name, sep = ""))){
    print("results directory already exists")
  }  else{
    dir.create(result_dir_name)
    print("results folder has been created")
  }
  
  # create plot directory
  if (dir.exists(paste(getwd(), plots_dir_name, sep = ""))){
    print("plots directory already exists")
  }  else{
    dir.create(plots_dir_name)
    print("plots folder has been created")
  }
  
  # read MS data from mzML file list
  msd <- readMSData(files = files,
                    pdata = new("NAnnotatedDataFrame", phenodata),
                    mode = "onDisk",
                    centroided = TRUE)
  
  # restrict data to length of chromatography run
  msd <- filterRt(msd, c(0, chrom_run_sec))
  
  # restrict data to MS1
  if (is.null(msLevel)) {
    msLevel <- 1
  } else {
    msLevel <- msLevel
  }
  
  msd <- filterMsLevel(msd, msLevel = msLevel)
  
  # create csv of raw data
  write.csv(fData(msd), file=paste(result_dir_name, filename = "_raw_data.csv", sep = ""), row.names=FALSE)
  
  save(msd, file = paste(result_dir_name, "/msd.RData", sep = ""))
  return(msd)
  
}

# ---- chromatograms ----
chromatogram_qc <- function(msd, phenodata, result_dir_name, plots_dir_name) {
  ### create quality control
  # create color palette
  color_pal <- brewer.pal(n = 9, name = "Set1")  
  
  # create color vector from phenodata
  uniq_sample_group <- data.frame(unique(phenodata$sample_group))
  uniq_sample_group$color <- color_pal[1:nrow(uniq_sample_group)]
  
  for (i in 1:nrow(uniq_sample_group)) {
    for (j in 1:nrow(phenodata)) {
      if (uniq_sample_group[i,1] %in% phenodata$sample_group[j]){
        phenodata$color[j] <- uniq_sample_group$color[i]
      }
    }
  }
  
  color <- phenodata$color

  # create chromatogram object for quality control plots
  chromas <- chromatogram(msd, aggregationFun="max", msLevel = 1)

  # Plot base peak chromatograms based on phenodata groups
  jpeg(filename = paste(plots_dir_name, "/chromas_bpc.jpeg", sep = ""), width = 2000, height = 1200, quality = 100, bg = "white")
  par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
  plot(chromas, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col = color)
  legend("topleft", bty="n", pt.cex=3, cex=1.5, y.intersp=0.7, text.width=0.5, pch=20,
         col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
  dev.off()

  # Plot total ion current chromatogram
  jpeg(filename = paste(plots_dir_name, "/chromas_tic.jpeg", sep = ""), width = 2000, height = 1200, quality = 100, bg = "white")
  par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0), cex.axis=1.5, cex=2, cex.lab=2, cex.main=2)
  tics <- split(tic(msd), f=fromFile(msd))
  boxplot(tics, col=color, ylab="intensity", xlab="sample", main="Total ion current", outline = FALSE)
  legend("topleft", bty="n", pt.cex=2, cex=1, y.intersp=0.7, text.width=0.5, pch=20,
         col= unique(color), legend= unique(msd@phenoData@data[["sample_group"]]))
  dev.off()

  # Qugrouping/binning for peak detection, based on similarity of the base peak chromatogram
  chromas_bin <- MSnbase::bin(chromas, binSize=2)
  chromas_bin_cor <- cor(log2(do.call(cbind, lapply(chromas, intensity)))) # transformation with log
  colnames(chromas_bin_cor) <- rownames(chromas_bin_cor) <- msd$sample_name
  chromas_bin_cor[is.na(chromas_bin_cor)] <- 0

  # representing the data in a heatmap for general overview
  jpeg(filename = paste(plots_dir_name, "/heatmap_chromas_bin.jpeg", sep = ""), width = 500, height = 500, quality = 100, bg = "white")
  par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(7,0,0,7), cex.axis=0.9, cex=0.6)
  heatmap(chromas_bin_cor)
  dev.off()
  
  # save
  save(chromas_bin, file = paste(result_dir_name, "/chromas_bin.RData", sep = ""))
  return(chromas_bin)

}


# ---- peak_detection -----
peak_detection <- function(msd, CentWaveParam = NULL, result_dir_name, plots_dir_name) {
  # msd = prepared ms object from data_preparation function
  # CentWaveParam = CentWaVeParam class for parameters for peak detection algorithm
  # 
  
  # define peak detection parameters, if default
  if (is.null(CentWaveParam)){
    CentWaveParam <- CentWaveParam(ppm=25, mzCenterFun="wMean", peakwidth=c(20, 50), 
                                    prefilter=c(3, 100), mzdiff=0.0155, snthresh=10, noise=0, 
                                    integrate=1L, firstBaselineCheck=TRUE, verboseColumns=FALSE, 
                                    fitgauss=FALSE, roiList=list(), roiScales=numeric())
  } else {
    CentWaveParam <- CentWaveParam
  }
  
  # perform peak detection
  ms1_data <- findChromPeaks(msd, param=CentWaveParam)
  
  # create summary file .csv
  ms1_summary <- lapply(split.data.frame(chromPeaks(ms1_data), 
                                                  f=chromPeaks(ms1_data)[, "sample"]), 
                                 FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
  ms1_summary <- do.call(rbind, ms1_summary)
  rownames(ms1_summary) <- basename(fileNames(ms1_data))
  
  write.csv(as.data.frame(table(msLevel(ms1_data))), file=paste(result_dir_name, "/ms1_data.csv", sep = ""), row.names=FALSE)
  
  # # create quality control plots
  # jpeg(filename = paste(plots_dir_name, "/ms1_data.jpeg"), width = 1000, height = 1000, quality = 100, bg = "white")
  # par(mfrow=c(1,1), mar=c(5,18,4,1), oma=c(0,0,0,0), cex.axis=1, cex=2, cex.lab=2, cex.main=2)
  # plotChromPeakImage(ms1_data, main="Frequency of identified peaks per RT", binSize = 20)
  # dev.off()
  # 
  # save detected object
  save(ms1_data, file = paste(result_dir_name, "/ms1_data.RData", sep = ""))
  return(ms1_data)
}

# ----grouping_1 ----
grouping_1 <- function(ms1_data, PeakDensityParam_gr = NULL, result_dir_name) {
  # ms1_data = object from peak detection step
  # PeakDensityParam = parameters for peak grouping
  
  # define parameters for grouping
  if (is.null(PeakDensityParam_gr)) {
    PeakDensityParam_gr <- PeakDensityParam(sampleGroups = ms1_data$sample_group, minFraction = 0.7, bw = 2.5)
  } else {
    PeakDensityParam_gr <- PeakDensityParam
  }
  
  # grouping
  ms1_data <- groupChromPeaks(ms1_data, param = PeakDensityParam_gr)
  
  # save object
  save(ms1_data, file = paste(result_dir_name, "/ms1_data.RData", sep = ""))
  return(ms1_data)
}

# ----rt_correction ----
rt_correction <- function(ms1_data, plots_dir_name, result_dir_name, PeakGroupsParam_rt, chromas_msd) {
  # ms1_data = 
  # PeakGroupsParam = parameters for retention time correction
  
  # parameters for retention time correction 
  if (is.null(PeakGroupsParam_rt)) {
    PeakGroupsParam_rt <- PeakGroupsParam(minFraction = 0.7, smooth="loess",span=0.5,family="gaussian")
  } else {
    PeakGroupsParam_rt <- PeakGroupsParam
  }
  
  # retention time correction
  ms1_data <- adjustRtime(ms1_data, param=PeakGroupsParam_rt)
  
  
  # Plot the difference of raw and adjusted retention times
  jpeg(filename = paste(plots_dir_name, "/ms1_raw_adjusted.jpeg", sep = ""), width = 1000, height = 2000, quality = 100, bg = "white")
  par(mfrow=c(2,1), mar=c(5,6,4,1), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
  plot(chromas_msd, peakType="none", main="Raw chromatograms")
  plotAdjustedRtime(ms1_data, lwd=2, main="Retention Time correction")
  par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(1,1,0,0), cex.axis=0.9, cex=0.8)
  dev.off()
  
  # save object
  save(ms1_data, file = paste(result_dir_name, "/ms1_data.RData", sep = ""))
  return(ms1_data)
}


# ----grouping_2----
grouping_2 <- function(ms1_data, PeakDensityParam_gr = NULL, result_dir_name) {
  # ms1_data = object from peak detection step
  # PeakDensityParam = parameters for peak grouping
  
  # define parameters for grouping
  if (is.null(PeakDensityParam_gr)) {
    PeakDensityParam_gr <- PeakDensityParam(sampleGroups = ms1_data$sample_group, minFraction = 0.7, bw = 2.5)
  } else {
    PeakDensityParam_gr <- PeakDensityParam_gr
  }
  
  # grouping
  ms1_data <- groupChromPeaks(ms1_data, param = PeakDensityParam_gr)
  
  # save object
  save(ms1_data, file = paste(result_dir_name, "/ms1_data_peak_det.RData", sep = ""))
  return(ms1_data)

}


# ----feature_extraction ----
feature_extraction <- function(ms1_data, result_dir_name) {
  #
  #
  
  # Build feature matrix
  ms1_matrix <- featureValues(ms1_data, method="medret", value="into")
  colnames(ms1_matrix) <- ms1_data@phenoData@data[["sample_name"]]
  write.csv(ms1_matrix, paste(result_dir_name, "/ms1_matrix.csv", sep = ""))
  
  # transpose feature table
  feature_list <- t(ms1_matrix)
  
  # Build feature summary
  ms1_summary <- featureSummary(ms1_data)
  save(ms1_summary, file = paste(result_dir_name, "/ms1_summary.RData", sep = ""))
  
  ms1_definition <- featureDefinitions(ms1_data)
  save(ms1_definition, file = paste(result_dir_name, "/ms1_definition.RData", sep = ""))
  
  # save feature table
  write.csv(feature_list, file = paste(result_dir_name, "/feature_list.csv", sep = ""))
  return(feature_list)
}


# ----feature_transformation----
feature_transformation <- function(feature_list, result_dir_name, plots_dir_name) {
  ### Transform data
  feature_list <- log2(feature_list)
  
  # change 0 to small value to distinguish between values of 1 and NA
  feature_list[which(feature_list == 0)] <- 0.01
  
  # Missing value imputation
  feature_list[which(is.na(feature_list))] <- 0
  
  # Plot histogram
  # jpeg(filename = paste(plots_dir_name, "/feat_list_hist.jpeg"), width = 1000, height = 1000, quality = 100, bg = "white")
  # hist(as.numeric(feature_list), main="Histogram of feature table")
  # dev.off()
  # 
  # save as csv
  write.csv(feature_list, file=paste(result_dir_name, "/feature_list.csv", sep = ""))
  return(feature_list)
  
}

# ----bina_list_creation-----
bina_list_creation <- function(ms1_matrix = NULL, intensity_cutoff, result_dir_name) {
  # ms1_matrix = csv created before, load from csv
  if (is.null(ms1_matrix)){
    ms1_matrix <- read.csv(file = paste(result_dir_name, "/ms1_matrix.csv", sep = ""))
    rownames(ms1_matrix) <- ms1_matrix[,1]
    ms1_matrix <- ms1_matrix[,-1]
  } else {
    ms1_matrix <- ms1_matrix
  }
  
  # Create single 0/1 matrix
  binary_list <- t(ms1_matrix)
  binary_list[is.na(binary_list)] <- 1
  binary_list <- log2(binary_list)
  binary_list[binary_list < intensity_cutoff] <- 0
  binary_list[binary_list != 0] <- 1
  
  # save as csv
  write.csv(binary_list, file=paste(result_dir_name, "/binary_list.csv", sep = ""))
  return(binary_list)
}

