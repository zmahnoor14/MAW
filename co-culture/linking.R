#### First link the feature table with inclusion list

##### separate inclusion lists ENDO
endo_inc <- read.csv("/KSS_210324_ENDO.csv", sep =";")
colnames(endo_inc) <- c("Mass",	"Formula",	"Formula_type",	"Species",	"CS",	"Polarity",	"Start",	"End",	"CE",	"Comment")
endo_inc["Start"] <- endo_inc["Start"]*60
endo_inc["End"] <- endo_inc["End"]*60
neg <- endo_inc[endo_inc["Polarity"]=="Negative" , ]
pos <- endo_inc[endo_inc["Polarity"]=="Positive" , ]
write.csv(neg, "KSS_210324_ENDO_neg.csv")
write.csv(pos, "KSS_210324_ENDO_pos.csv")

##### separate inclusion lists EXO
exo_inc <- read.csv("/KSS_210324_EXO.csv", sep =";")
colnames(exo_inc) <- c("Mass",	"Formula",	"Formula_type",	"Species",	"CS",	"Polarity",	"Start",	"End",	"CE",	"Comment")
exo_inc["Start"] <- exo_inc["Start"]*60
exo_inc["End"] <- exo_inc["End"]*60
neg <- exo_inc[exo_inc["Polarity"]=="Negative" , ]
pos <- exo_inc[exo_inc["Polarity"]=="Positive" , ]
write.csv(neg, "KSS_210324_EXO_neg.csv")
write.csv(pos, "KSS_210324_EXO_pos.csv")

#---------------------------------------------------------------------------------------------------------------------------

##### Linking #####


linking_with_inc_list <- function(ft_csv, inc_csv, cond, mass_col, rtmin_col, rtmax_col){
  # make feature table smaller by excluding features not in MS2 inclusion list
  ft_tbl <- read.csv(ft_csv)
  inc_tbl <- read.csv(inc_csv)

  inc_tbl['feat_inc'] <- NA
  inc_tbl["Comment"]<- NA
  
  for (i in 1:nrow(inc_tbl)){

    for (j in 1:nrow(ft_tbl)){
      if (inc_tbl[i, mass_col]>ft_tbl[j, "mzmin"] && ft_tbl[j, "mzmax"]>inc_tbl[i, mass_col]){
        if (ft_tbl[j, "rtmed"]>inc_tbl[i, rtmin_col] && inc_tbl[i, rtmax_col]>ft_tbl[j, "rtmed"]){
          inc_tbl[i, 'feat_inc'] <- ft_tbl[j, "feat_id"]
        }
      }
      else if(round(inc_tbl[i, mass_col],4)>=round(ft_tbl[j, "mzmin"],4) && round(ft_tbl[j, "mzmax"],4)>=round(inc_tbl[i, mass_col],4)){
        if (ft_tbl[j, "rtmed"]>inc_tbl[i, rtmin_col] && inc_tbl[i, rtmax_col]>ft_tbl[j, "rtmed"]){
          inc_tbl[i, 'feat_inc'] <- ft_tbl[j, "feat_id"]
        }
      }
      else if(round(inc_tbl[i, mass_col],3)>=round(ft_tbl[j, "mzmin"],3) && round(ft_tbl[j, "mzmax"],3)>=round(inc_tbl[i, mass_col],3)){
        if (ft_tbl[j, "rtmed"]>inc_tbl[i, rtmin_col] && inc_tbl[i, rtmax_col]>ft_tbl[j, "rtmed"]){
          inc_tbl[i, 'feat_inc'] <- ft_tbl[j, "feat_id"]
        }
      }
    }
  }
  
  
  for (i in 1:nrow(inc_tbl)){
    if (is.na(inc_tbl[i, "feat_inc"])){
      ft_tbl_ids <- c()
      ft_rt_med <- c()
      for (j in 1:nrow(ft_tbl)){
        if(round(inc_tbl[i, mass_col],2)>=round(ft_tbl[j, "mzmin"],2) && round(ft_tbl[j, "mzmax"],2)>=round(inc_tbl[i, mass_col],2)){
          ft_tbl_ids <- c(ft_tbl_ids, ft_tbl[j, "feat_id"])
          ft_rt_med <- c(ft_rt_med, ft_tbl[j, "rtmed"])
        }
      }
      ft_tbl_rt_sh <- cbind(ft_tbl_ids, ft_rt_med)
      print(nrow(ft_tbl_rt_sh))
      if (!(is.null(ft_tbl_rt_sh))){
        if (nrow(ft_tbl_rt_sh) == 1){
          inc_tbl[i, 'feat_inc'] <- ft_tbl_ids
          inc_tbl[i, "Comment"] <- paste("The ONLY rtmed for this premz is:", ft_tbl_rt_sh[1,"ft_rt_med"])
        }
        else if (nrow(ft_tbl_rt_sh)>1){
          # Calculate the absolute differences between each number and the range boundaries
          differences <- abs(as.numeric(ft_tbl_rt_sh[,"ft_rt_med"]) - inc_tbl[i, rtmin_col])
          differences <- pmin(differences, abs(as.numeric(ft_tbl_rt_sh[,"ft_rt_med"]) - inc_tbl[i, rtmax_col]))
          
          # Find the index of the number with the minimum difference
          closest_index <- which.min(differences)
          inc_tbl[i, 'feat_inc'] <- ft_tbl_rt_sh[closest_index, "ft_tbl_ids"]
          # Get the closest number
          closest_number <- as.numeric(ft_tbl_rt_sh[,"ft_rt_med"])[closest_index]
          
          # Print the result
          inc_tbl[i, "Comment"] <- paste("The rtmed closest to the range is:", closest_number)
        }
      }
    }
  }
  
  write.csv(inc_tbl,paste("ms_tbl_", cond, ".csv", sep = ""))
  return(inc_tbl)
}

link_ann_inc <- function(inc_tbl, ann_tbl, cond){
  ann_tbl <- read.csv(ann_tbl) # results from MAW
  inc_tbl <- read.csv(inc_tbl) # result from previous function
  ann_tbl['feat_inc'] <- NA
  ann_tbl["Comment"]<- NA
  for (i in 1:nrow(ann_tbl)) {
    for (j in 1:nrow(inc_tbl)) {
      if (round(ann_tbl[i, "premz"], 4) == round(inc_tbl[j, "Mass"], 4)){
        if (ann_tbl[i, "rtmed"]>inc_tbl[j, "Start"] && inc_tbl[j,"End"]>ann_tbl[i, "rtmed"]){
          x <- x + 1
          ann_tbl[i, "Comment"] <- inc_tbl[j, "Comment"]
          ann_tbl[i, "feat_inc"] <- inc_tbl[j, "feat_inc"]
        }
      }
      else if (round(ann_tbl[i, "premz"], 3) == round(inc_tbl[j, "Mass"], 3)){
        if (ann_tbl[i, "rtmed"]>inc_tbl[j, "Start"] && inc_tbl[j,"End"]>ann_tbl[i, "rtmed"]){
          x <- x + 1
          ann_tbl[i, "Comment"] <- inc_tbl[j, "Comment"]
          ann_tbl[i, "feat_inc"] <- inc_tbl[j, "feat_inc"]
        }
      }
    }
  }

  for (i in 1:nrow(ann_tbl)){
    if (!(grepl("FT", ann_tbl[i, "feat_inc"]))){
      ft_tbl_ids <- c()
      ft_rt_min <- c()
      ft_rt_max <- c()
      for (j in 1:nrow(inc_tbl)){
        if(round(ann_tbl[i, "premz"], 3) == round(inc_tbl[j, "Mass"], 3)){
          ft_tbl_ids <- c(ft_tbl_ids, inc_tbl[j, "feat_inc"])
          ft_rt_min <- c(ft_rt_min, inc_tbl[j, "Start"])
          ft_rt_max <- c(ft_rt_max, inc_tbl[j, "End"])
        }
      }
      ft_tbl_rt_sh <- cbind(ft_tbl_ids, ft_rt_min, ft_rt_max)
      print(nrow(ft_tbl_rt_sh))
      if (nrow(ft_tbl_rt_sh) == 1){
        ann_tbl[i, 'feat_inc'] <- ft_tbl_ids
        ann_tbl[i, "Comment"] <- paste("The ONLY range for this ann_premz is:", ft_tbl_rt_sh[1,"ft_rt_min"], "-", ft_tbl_rt_sh[1,"ft_rt_max"], sep = "")
      }
      else if (nrow(ft_tbl_rt_sh)>1){
        # Calculate the absolute differences between each number and the range boundaries
        differences <- abs(as.numeric(ann_tbl[i,"rtmed"]) - as.numeric(ft_tbl_rt_sh[, "ft_rt_min"]))
        differences <- pmin(differences, abs(as.numeric(ann_tbl[i,"rtmed"]) - as.numeric(ft_tbl_rt_sh[, "ft_rt_max"])))
        # Find the index of the number with the minimum difference
        closest_index <- which.min(differences)
        ann_tbl[i, 'feat_inc'] <- ft_tbl_rt_sh[closest_index, "ft_tbl_ids"]
        # Get the closest number
        closest_min <- as.numeric(ft_tbl_rt_sh[,"ft_rt_min"])[closest_index]
        closest_max <- as.numeric(ft_tbl_rt_sh[,"ft_rt_max"])[closest_index]
        # Print the result
        ann_tbl[i, "Comment"] <- paste("The rtmed closest to range is:", closest_min, "-", closest_min)
      }
    }
  }


  write.csv(ann_tbl, paste(cond, "_ann.csv", sep = ""))
  return(ann_tbl)
}

assign_sig_feat <- function(sig_feat, ann_tbl_modified, cond, pol, feat_col, species){
  
  sig_feat <- read.csv(sig_feat)
  ann_tbl <- read.csv(ann_tbl_modified)
  sig_feat_pos <- sig_feat[grepl(pol, sig_feat[,feat_col]), ]
  
  sig_feat_pos["sig_feat_ids"] <-gsub(pol, '', sig_feat_pos[,feat_col])
  sig_feat_pos["sig_feat_ids"]
  ann_tbl[paste(species, "Co", sep = "")] <- NA
  ann_tbl[paste(species, "Mono", sep = "")] <- NA

  for (i in 1:nrow(ann_tbl)){
    if (!(is.na(as.character(ann_tbl[i, "feat_inc"])))){
      for (j in 1:nrow(sig_feat_pos)){
        if (as.character(ann_tbl[i, "feat_inc"]) == as.character(sig_feat_pos[j, "sig_feat_ids"])){
          if ("Co" == as.character(sig_feat_pos[j, "Condition"])){
            ann_tbl[i, paste(species, "Co", sep = "")] <- "present"
          }
          else if ("Mono" == as.character(sig_feat_pos[j, "Condition"])){
            ann_tbl[i, paste(species, "Mono", sep = "")] <- "present" 
          }
        }
      }
    }
  }
  write.csv(ann_tbl, paste(cond, "_mergedResults-with-one-Candidates_sig_feat_for_only_inclusion.csv", sep = ""))
  return(ann_tbl)
}



#---------------------------------------------------------------------------------------------------------------------------

##### script for exo_pos #####



setwd("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking")

ft_csv = "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/feature_info_EXO_pos.csv"
inc_csv = "KSS_210324_EXO_pos.csv"
cond = "exo_pos"
mass_col = "Mass"
rtmin_col = "Start"
rtmax_col = "End"
ann_tbl <- "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/exo_pos_sirius_mergedResults-with-one-Candidates.csv"
ann_tbl_modified <- paste(cond, "_ann.csv", sep ="")
sig_feat <- "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/new_sel_pls_features_Pp_bina.csv"
feat_col <- "Feature"
pol = "_pos" # or "_neg"
species = "Pp"

inc_tbl1 <- linking_with_inc_list(ft_csv, inc_csv, cond, mass_col, rtmin_col, rtmax_col)
inc_tbl1 

ann_tbl1 <- link_ann_inc(inc_tbl = paste("ms_tbl_", cond, ".csv", sep = ""), ann_tbl, cond)
ann_tbl1

sig_features <- assign_sig_feat(sig_feat, 
                            ann_tbl_modified, 
                            cond, 
                            pol, 
                            feat_col,
                            species)

no_of_sig_feats_smco<- nrow(sig_features[!(is.na(sig_features["SmCo"])), ])
no_of_sig_feats_smco
no_of_sig_feats_ppco<- nrow(sig_features[!(is.na(sig_features["PpCo"])), ])
no_of_sig_feats_ppco
no_of_sig_feats_smmono<- nrow(sig_features[!(is.na(sig_features["SmMono"])), ])
no_of_sig_feats_smmono
no_of_sig_feats_ppmono<- nrow(sig_features[!(is.na(sig_features["PpMono"])), ])
no_of_sig_feats_ppmono



sgg <- sig_features[!(is.na(sig_features["origin"])), ]
sgg
##### script for exo_neg #####

setwd("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking")

ft_csv = "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/feature_info_EXO_neg.csv"
inc_csv = "KSS_210324_EXO_neg.csv"
cond = "exo_neg"
mass_col = "Mass"
rtmin_col = "Start"
rtmax_col = "End"
ann_tbl <- "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/exo_neg_sirius_mergedResults-with-one-Candidates.csv"
ann_tbl_modified <- paste(cond, "_ann.csv", sep ="")
sig_feat <- "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/new_sel_pls_features_Sm_bina.csv"
feat_col <- "Feature"
pol = "_neg"
species = "Sm"

inc_tbl1 <- linking_with_inc_list(ft_csv, inc_csv, cond, mass_col, rtmin_col, rtmax_col)
inc_tbl1 

ann_tbl1 <- link_ann_inc(inc_tbl = paste("ms_tbl_", cond, ".csv", sep = ""), ann_tbl, cond)
ann_tbl1

sig_features <- assign_sig_feat(sig_feat, 
                                ann_tbl_modified, 
                                cond, 
                                pol, 
                                feat_col,
                                species)

no_of_sig_feats_smco<- nrow(sig_features[!(is.na(sig_features["SmCo"])), ])
no_of_sig_feats_smco
no_of_sig_feats_ppco<- nrow(sig_features[!(is.na(sig_features["PpCo"])), ])
no_of_sig_feats_ppco
no_of_sig_feats_smmono<- nrow(sig_features[!(is.na(sig_features["SmMono"])), ])
no_of_sig_feats_smmono
no_of_sig_feats_ppmono<- nrow(sig_features[!(is.na(sig_features["PpMono"])), ])
no_of_sig_feats_ppmono


##### script for endo_neg #####



setwd("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking")

ft_csv = "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/feature_info_ENDO_neg.csv"
inc_csv = "KSS_210324_ENDO_neg.csv"
cond = "endo_neg"
mass_col = "Mass"
rtmin_col = "Start"
rtmax_col = "End"
ann_tbl <- "endo_neg_mergedResults-with-one-Candidates.csv"
ann_tbl_modified <- paste(cond, "_ann.csv", sep ="")
sig_feat <- "sel_pls_feat_ENDO.csv"
feat_col <- "sel_pls_features_ENDO"
origin_col <- "origin"
pol = "_neg"

inc_tbl3 <- linking_with_inc_list(ft_csv, inc_csv, cond, mass_col, rtmin_col, rtmax_col)
inc_tbl3 

ann_tbl3 <- link_ann_inc(inc_tbl = paste("ms_tbl_", cond, ".csv", sep = ""), ann_tbl, cond)
ann_tbl3

sig_features3 <- assign_sig_feat(sig_feat, 
                                ann_tbl_modified, 
                                cond, 
                                pol, 
                                feat_col, 
                                origin_col)

no_of_sig_feats3<- nrow(sig_features3[!(is.na(sig_features3["origin"])), ])
no_of_sig_feats3

##### script for endo_pos #####



setwd("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking")

ft_csv = "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/feature_info_ENDO_pos.csv"
inc_csv = "KSS_210324_ENDO_pos.csv"
cond = "endo_pos"
mass_col = "Mass"
rtmin_col = "Start"
rtmax_col = "End"
ann_tbl <- "endo_pos_mergedResults-with-one-Candidates.csv"
ann_tbl_modified <- paste(cond, "_ann.csv", sep ="")
sig_feat <- "sel_pls_feat_ENDO.csv"
feat_col <- "sel_pls_features_ENDO"
origin_col <- "origin"
pol = "_pos"

inc_tbl4 <- linking_with_inc_list(ft_csv, inc_csv, cond, mass_col, rtmin_col, rtmax_col)
inc_tbl4 

ann_tbl4 <- link_ann_inc(inc_tbl = paste("ms_tbl_", cond, ".csv", sep = ""), ann_tbl, cond)
ann_tbl4

sig_features4 <- assign_sig_feat(sig_feat, 
                                 ann_tbl_modified, 
                                 cond, 
                                 pol, 
                                 feat_col, 
                                 origin_col)

no_of_sig_feats4 <- nrow(sig_features4[!(is.na(sig_features4["origin"])), ])
sig <- sig_features4[!(is.na(sig_features4["origin"])), ]


##### Linking conditions to the feature ids for endo neg#####

load("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/CoCulture/endo_neg_Results/MS1_ENDO_neg_peak_detection.RData")

ft_id <- ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@rownames
PpCoCu <- ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["CoCuPp"]]
PpMono <- ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["Pp"]]
SmCoCu <- ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["CoCuSm"]]
SmMono <- ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["Sm"]]
MB <- ms1_data_ENDO_neg@msFeatureData[["featureDefinitions"]]@listData[["MB"]]

int_feat_df <- cbind(ft_id, SmMono, SmCoCu, PpMono, PpCoCu, MB)

ann_feat_df <- read.csv("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/endo_neg_mergedResults-with-one-Candidates_sig_feat_for_only_inclusion.csv")

ann_feat_df["PpCoCu"] <- NA 
ann_feat_df["SmCoCu"] <- NA 
ann_feat_df["PpMono"] <- NA 
ann_feat_df["SmMono"] <- NA 
ann_feat_df["MB"] <- NA 

for (i in 1:nrow(ann_feat_df)){
  if (!(is.na(as.character(ann_feat_df[i, "feat_inc"])))){
    for (j in 1:nrow(int_feat_df)){
      if (as.character(ann_feat_df[i, "feat_inc"]) == as.character(int_feat_df[j, "ft_id"])){
        ann_feat_df[i, "PpCoCu"] <- int_feat_df[j, "PpCoCu"] 
        ann_feat_df[i, "SmCoCu"] <- int_feat_df[j, "SmCoCu"] 
        ann_feat_df[i, "PpMono"] <- int_feat_df[j, "PpMono"] 
        ann_feat_df[i, "SmMono"] <- int_feat_df[j, "SmMono"] 
        ann_feat_df[i, "MB"] <- int_feat_df[j, "MB"] 
      }
    }
  }
}
write.csv(ann_feat_df, "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/endo_neg_ann_origin_feat_df.csv")

##### Linking conditions to the feature ids for endo pos#####

load("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/CoCulture/endo_pos_Results/MS1_ENDO_pos_peak_detection.RData")

ft_id <- ms1_data_ENDO_pos@msFeatureData[["featureDefinitions"]]@rownames
PpCoCu <- ms1_data_ENDO_pos@msFeatureData[["featureDefinitions"]]@listData[["CoCuPp"]]
PpMono <- ms1_data_ENDO_pos@msFeatureData[["featureDefinitions"]]@listData[["Pp"]]
SmCoCu <- ms1_data_ENDO_pos@msFeatureData[["featureDefinitions"]]@listData[["CoCuSm"]]
SmMono <- ms1_data_ENDO_pos@msFeatureData[["featureDefinitions"]]@listData[["Sm"]]
MB <- ms1_data_ENDO_pos@msFeatureData[["featureDefinitions"]]@listData[["MB"]]

int_feat_df <- cbind(ft_id, SmMono, SmCoCu, PpMono, PpCoCu, MB)

ann_feat_df <- read.csv("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/endo_pos_mergedResults-with-one-Candidates_sig_feat_for_only_inclusion.csv")

ann_feat_df["PpCoCu"] <- NA 
ann_feat_df["SmCoCu"] <- NA 
ann_feat_df["PpMono"] <- NA 
ann_feat_df["SmMono"] <- NA 
ann_feat_df["MB"] <- NA 

for (i in 1:nrow(ann_feat_df)){
  if (!(is.na(as.character(ann_feat_df[i, "feat_inc"])))){
    for (j in 1:nrow(int_feat_df)){
      if (as.character(ann_feat_df[i, "feat_inc"]) == as.character(int_feat_df[j, "ft_id"])){
        ann_feat_df[i, "PpCoCu"] <- int_feat_df[j, "PpCoCu"] 
        ann_feat_df[i, "SmCoCu"] <- int_feat_df[j, "SmCoCu"] 
        ann_feat_df[i, "PpMono"] <- int_feat_df[j, "PpMono"] 
        ann_feat_df[i, "SmMono"] <- int_feat_df[j, "SmMono"] 
        ann_feat_df[i, "MB"] <- int_feat_df[j, "MB"] 
      }
    }
  }
}
write.csv(ann_feat_df, "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/endo_pos_ann_origin_feat_df.csv")

##### Linking conditions to the feature ids for exo neg#####

load("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/CoCulture/exo_neg_Results/MS1_EXO_neg_peak_detection.RData")

ft_id <- MS1_EXO_neg_peak_detection@msFeatureData[["featureDefinitions"]]@rownames
PpCoCu <- MS1_EXO_neg_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["CoCuPp"]]
PpMono <- MS1_EXO_neg_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["Pp"]]
SmCoCu <- MS1_EXO_neg_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["CoCuSm"]]
SmMono <- MS1_EXO_neg_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["Sm"]]
MB <- MS1_EXO_neg_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["MB"]]

int_feat_df <- cbind(ft_id, SmMono, SmCoCu, PpMono, PpCoCu, MB)

ann_feat_df <- read.csv("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/exo_neg_mergedResults-with-one-Candidates_sig_feat_for_only_inclusion.csv")

ann_feat_df["PpCoCu"] <- NA 
ann_feat_df["SmCoCu"] <- NA 
ann_feat_df["PpMono"] <- NA 
ann_feat_df["SmMono"] <- NA 
ann_feat_df["MB"] <- NA 

for (i in 1:nrow(ann_feat_df)){
  if (!(is.na(as.character(ann_feat_df[i, "feat_inc"])))){
    for (j in 1:nrow(int_feat_df)){
      if (as.character(ann_feat_df[i, "feat_inc"]) == as.character(int_feat_df[j, "ft_id"])){
        ann_feat_df[i, "PpCoCu"] <- int_feat_df[j, "PpCoCu"] 
        ann_feat_df[i, "SmCoCu"] <- int_feat_df[j, "SmCoCu"] 
        ann_feat_df[i, "PpMono"] <- int_feat_df[j, "PpMono"] 
        ann_feat_df[i, "SmMono"] <- int_feat_df[j, "SmMono"] 
        ann_feat_df[i, "MB"] <- int_feat_df[j, "MB"] 
      }
    }
  }
}
ann_feat_df
write.csv(ann_feat_df, "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/exo_neg_ann_origin_feat_df.csv")

##### Linking conditions to the feature ids for exo neg#####

load("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/CoCulture/exo_pos_Results/MS1_EXO_pos_peak_detection.RData")

ft_id2 <- MS1_EXO_pos_peak_detection@msFeatureData[["featureDefinitions"]]@rownames
PpCoCu2 <- MS1_EXO_pos_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["CoCuPp"]]
PpMono2 <- MS1_EXO_pos_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["Pp"]]
SmCoCu2 <- MS1_EXO_pos_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["CoCuSm"]]
SmMono2 <- MS1_EXO_pos_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["Sm"]]
MB2 <- MS1_EXO_pos_peak_detection@msFeatureData[["featureDefinitions"]]@listData[["MB"]]

int_feat_df2 <- cbind(ft_id, SmMono, SmCoCu, PpMono, PpCoCu, MB)

ann_feat_df2 <- read.csv("/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/exo_pos_mergedResults-with-one-Candidates_sig_feat_for_only_inclusion.csv")

ann_feat_df2["PpCoCu"] <- NA 
ann_feat_df2["SmCoCu"] <- NA 
ann_feat_df2["PpMono"] <- NA 
ann_feat_df2["SmMono"] <- NA 
ann_feat_df2["MB"] <- NA 

for (i in 1:nrow(ann_feat_df2)){
  if (!(is.na(as.character(ann_feat_df2[i, "feat_inc"])))){
    for (j in 1:nrow(int_feat_df2)){
      if (as.character(ann_feat_df2[i, "feat_inc"]) == as.character(int_feat_df2[j, "ft_id"])){
        ann_feat_df2[i, "PpCoCu"] <- int_feat_df2[j, "PpCoCu"] 
        ann_feat_df2[i, "SmCoCu"] <- int_feat_df2[j, "SmCoCu"] 
        ann_feat_df2[i, "PpMono"] <- int_feat_df2[j, "PpMono"] 
        ann_feat_df2[i, "SmMono"] <- int_feat_df2[j, "SmMono"] 
        ann_feat_df2[i, "MB"] <- int_feat_df2[j, "MB"] 
      }
    }
  }
}
write.csv(ann_feat_df2,  "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/exo_pos_ann_origin_feat_df.csv")



