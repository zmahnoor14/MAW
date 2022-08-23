#! /usr/bin/Rscript

#' @title SIRIUS Adducts and Results Post-processing
#'
#' @description
#'
#' This function performs first part of SIRIUS Results post processing
#' i.e. finding the top candidate among a list of candidates; extracts top
#' 5 candidates SMILES; to calculate maximum common substrcuture in postprocessing part2

#' @param x is the result directory for the .mzml file
#'
#' @param SL either TRUE or FALSE (depends is the user has used 
#' inhouse library or here called a suspect list)

#' @return
#' 
#' a csv file containing features and SIRIUS ANNOTATION results
#' named MS1DATAsirius.csv
#' 
#' @author Mahnoor Zulfiqar
#' 
#' @examples
#' 
#' sirius_param(x = /usr/project/file1/',
#'                      SL = TRUE,


# ---------- Preparations ----------
# Load libraries
library(stringr)
library(readr)
library(dplyr)

# ---------- Arguments and user variables ----------
args <- commandArgs(trailingOnly=TRUE)
#print(args)

x <- as.character(args[1])

SL <- as.logical(args[2])

# ---------- sirius_postprocess ----------

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
            # file path for formula candidate from Param json folder
            for_canS <- paste(list.dirs(parameter_json[i,'SL_Param'])[2], '/formula_candidates.tsv', sep = '')
            # file path for strcuture candidate from Suspect List json folder
            str_can <- paste(list.dirs(parameter_json[i,'Param'])[2], '/structure_candidates.tsv', sep = '')
            # file path for formula candidate from Param json folder
            for_can <- paste(list.dirs(parameter_json[i,'Param'])[2], '/formula_candidates.tsv', sep = '')
            
            # if the strcuture candidate file exists
            if (file.exists(str_canS) && file.exists(str_can)){
            
                # read the corresponding structure and formula candidate files
                str_canSL <- as.data.frame(read_tsv(str_canS))
                for_canSL <- as.data.frame(read_tsv(for_canS))
                
                # read the corresponding structure and formula candidate files
                str_canP <- as.data.frame(read_tsv(str_can))
                for_canP <- as.data.frame(read_tsv(for_can))
            
                # if the strcuture candidate file contains 1 or more rows, it has detected a candidate from suspect list, add relevant info
                
                if (nrow(str_canSL) >= 1){
                     
                    if (str_canSL[1, 'CSI:FingerIDScore'] > str_canP[1, 'CSI:FingerIDScore']){
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
                    else{
                        
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
                # if it's empty, move onto the All DB result folder called PARAM here
                else{
                    
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



sirius_postprocess(x, SL)

