### This statistics file shows the exact analyses that were performed for the data at hand in the co-culture experiment

# Statistics analysis combining the feature tables of ENDO pos and ENDO neg
# Data from MS1_ENDO_neg and MS1_ENDO_pos scripts

###---- library ----
# Load libraries
#library(parallel)               # Detect number of cpu cores
#library(foreach)                # For multicore parallel
#library(doMC)                   # For multicore parallel
library(RColorBrewer)           # For colors
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(vegan)                  # For Statistics/Varpart
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(circlize)               # For sunburst plot
library(plotrix)                # For sunburst plot
library(caret)                  # Swiss-army knife for statistics
library(pROC)                   # Evaluation metrics
library(PRROC)                  # Evaluation metrics
library(multiROC)               # Evaluation metrics
library(chemodiv)               # Chemodiversity (Petren 2022)
#library(rcdk)                   # CDK
#library(rinchi)                 # Converting SMILES to InchiKey
library(plotly)                 # For creating html plots
library(htmlwidgets)            # For creating html plots
#library(shiny)                  # HTML in R
#library(sunburstR)              # HTML-sunburst plots
library(heatmaply)              # HTML heatmaps
library(stringr)              
library(randomForest)           # Random forest
library(MetaboAnalystR)         # Random forest
#library(iESTIMATE)
source("https://raw.githubusercontent.com/ipb-halle/iESTIMATE/main/R/_functions.r")

########## Import objects ########

# import feature tables
feature_table_pos <- read.csv("endo_pos/endo_pos_Results/feature_list_ENDO_pos.csv", row.names = 1)
# load(file = "ENDO_stats_Results/MS_ENDO_pos_peak_detection.RData")
feature_table_neg <- read.csv("endo_neg/endo_neg_Results/feature_list_ENDO_neg.csv", row.names = 1)
# load(file = "ENDO_neg_Results/MS_ENDO_pos_peak_detection.RData")

# name the features per polarity
colnames(feature_table_pos) <- paste(colnames(feature_table_pos), "pos", sep = "_")
colnames(feature_table_neg) <- paste(colnames(feature_table_neg), "neg", sep = "_")

# combine the feature tables of neg and pos into one
feat_list_ENDO <- cbind(feature_table_pos, feature_table_neg)

# import binary tables
binary_table_pos <- read.csv("endo_pos/endo_pos_Results/bina_list_ENDO_pos.csv", row.names = 1)
binary_table_neg <- read.csv("endo_neg/endo_neg_Results/bina_list_ENDO_neg.csv", row.names = 1)

# name the features per polarity
colnames(binary_table_pos) <- paste(colnames(binary_table_pos), "pos", sep = "_")
colnames(binary_table_neg) <- paste(colnames(binary_table_neg), "neg", sep = "_")

# combine the feature tables of neg and pos into one
bina_list_ENDO <- cbind(binary_table_pos, binary_table_neg)

###----Creation of Phenodata/Metadata----
# create vector with sample classes according to culture information sheet
samp_groups <- c("CoCuPp", "CoCuSm", "CoCuSm", "CoCuPp", "CoCuPp", "CoCuSm",
                 rep(x = "Sm", times = 8), 
                 rep(x = "Pp", times = 8),
                 "CoCuSm", "CoCuPp", "MB")

CoCuPp_des <- "Co-culture sample from Prymnesium parvum"
CoCuSm_des <- "Co-culture sample from Skeletonema marinoi"
Sm_des <- "Mono-culture sample from Skeletonema marinoi"
Pp_des <- "Mono-culture sample from Prymnesium parvum"
MB_des <- "Media Blank"

# create vector with sample classes according to culture information sheet
samp_groups_description <- c(CoCuPp_des, CoCuSm_des, CoCuSm_des, CoCuPp_des, CoCuPp_des, CoCuSm_des,
                             rep(x = Sm_des, times = 8), 
                             rep(x = Pp_des, times = 8),
                             CoCuSm_des, CoCuPp_des, MB_des)


# create vector with colors 
CoCuPp1 <- ("#3E134F")
CoCuSm1 <- rep("#F36E35", 2)
CoCuPp2 <- rep("#3E134F", 2)
CoCuSm2 <- ("#F36E35")
Sm <- rep("#F8B83C", 8)
Pp <- rep("#C53270", 8)
CoCuSm3 <- ("#F36E35")
CoCuPp3 <- ("#3E134F")
MB <- rep("#040404", 1)

color <- c(CoCuPp1, CoCuSm1, CoCuPp2, CoCuSm2, Sm, Pp, CoCuSm3, CoCuPp3, MB)

# file names
MS1_ENDO_names <- rownames(feat_list_ENDO)

# create phenodata based on culture type
pheno_data_ENDO <- data.frame(sample_name = MS1_ENDO_names, sample_group = samp_groups, samp_groups_description = samp_groups_description)
pheno_color_ENDO <- data.frame(color)

# set rownames of feature table
for (i in 1:24) {
  rownames(feat_list_ENDO)[i] <- gsub("coculture_pos", pheno_data_ENDO$sample_group[i], rownames(feat_list_ENDO)[i])
}

# set rownames of bina list
for (i in 1:24) {
  rownames(bina_list_ENDO)[i] <- gsub("coculture_pos", pheno_data_ENDO$sample_group[i], rownames(bina_list_ENDO)[i])
}


###----Preparations----
# create plot directory
if (dir.exists(paste(getwd(), "/ENDO_stats_plots/", sep = ""))){
  print("plots directory already exists")
}  else{
  dir.create("ENDO_stats_plots")
  print("plots folder has been created")
}

# create Results directory
if (dir.exists(paste(getwd(), "/ENDO_stats_Results/", sep = ""))){
  print("Results directory already exists")
}  else{
  dir.create("ENDO_stats_Results")
  print("Results folder has been created")
}



# define variables by which to analyse
mzml_pheno_samples <- samp_groups_description[1:24]
mzml_pheno_colors <- color[1:24]
# By Culture type and Species
mzml_pheno_samples_type <- samp_groups[1:24]
# By Culture type
mzml_pheno_origin_samples <- as.factor(origin_samp_groups <- c("CoCu", "CoCu", "CoCu", "CoCu", "CoCu", "CoCu",
                                                                   rep(x = "MonoCu", times = 8), rep(x = "MonoCu", times = 8),
                                                                   "CoCu", "CoCu"))
# By Species
mzml_pheno_origin_species <- as.factor(species_samp_groups <- c("P. parvum", "S. marinoi", "S. marinoi", "P. parvum", "P. parvum", "S. marinoi",
                                                                    rep(x = "S. marinoi", times = 8), rep(x = "P. parvum", times = 8),
                                                                    "S. marinoi", "P. parvum"))
# For plot legends
mzml_pheno_legend <-  c("Co-culture P.parvum", "Co-culture S.marinoi", "Co-culture S.marinoi", "Co-culture P.parvum", "Co-culture P.parvum", "Co-culture S.marinoi",
                        rep(x = "Mono-culture S.marinoi", times = 8), rep(x = "Mono-culture P.parvum", times = 8),
                        "Co-culture S.marinoi", "Co-culture P.parvum")


# overview dataframe
mzml_pheno <- data.frame(mzml_pheno_samples_type, mzml_pheno_origin_samples, mzml_pheno_origin_species, mzml_pheno_legend)


# remove media blank sample from analysis
feat_list_ENDO <- feat_list_ENDO[-25,]
bina_list_ENDO <- bina_list_ENDO[-25,]

# create objets for Statistical analysis
# comp_list 
comp_list <- feat_list_ENDO
comp_list <- comp_list[] # removes missing values
rownames(comp_list) <- paste(pheno_data_ENDO[1:24,1], pheno_data_ENDO[1:24,2], sep = "_" )
#colnames(comp_list) <- paste0(rownames(ms1_def_ENDO_pos)[which(ms1_def_ENDO_pos$has_ms2==1)], "_pos")

# bina_list
bina_list <- bina_list_ENDO
rownames(bina_list) <- paste(pheno_data_ENDO[1:24,1], pheno_data_ENDO[1:24,2], sep = "_" )
#colnames(bina_list) <- paste0(rownames(ms1_def_ENDO_pos)[which(ms1_def_ENDO_pos$has_ms2==1)], "_pos")

## save feat_list and bina_list
# save as csv
write.csv(feat_list_ENDO, file = "ENDO_stats_Results/feat_list_ENDO.csv")
write.csv(bina_list_ENDO, file = "ENDO_stats_Results/bina_list_ENDO.csv")

# ##################### create objects for analysis on species level ###########
# index for the species types
index_SM <- grep("Sm", pheno_data_ENDO$sample_group)
index_PP <- grep("Pp", pheno_data_ENDO$sample_group)

# divide the feature table by species type
feat_list_ENDO_Sm <- feat_list_ENDO[index_SM,]
feat_list_ENDO_Pp <- feat_list_ENDO[index_PP,]

# create objets for Statistical analysis
# comp_list species
comp_list_Sm <- feat_list_ENDO_Sm
comp_list_Sm <- comp_list_Sm[] # removes missing values

comp_list_Pp <- feat_list_ENDO_Pp
comp_list_Pp <- comp_list_Pp

# divide bina list by species type
bina_list_Sm <- bina_list_ENDO[index_SM,]
bina_list_Pp <- bina_list_ENDO[index_PP,]

# ############################## MS1 statistics  ##############################


# ---------- Diversity  ----------

# Create data frame
model_div             <- data.frame(features=apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div$richness    <- apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } )
model_div$menhinick   <- apply(X=bina_list, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div$shannon     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div$pielou      <- apply(X=scale(comp_list, center=FALSE, scale=FALSE), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list, species), index="chao")
model_div$simpson     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div$inverse     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div$fisher      <- apply(X=comp_list, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
#model_div$unique      <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) })
#model_div$hillfunc    <- as.numeric(unlist(calcDiv(comp_list, compDisMat=scales::rescale(as.matrix(dist(t(comp_list)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))
#rownames(model_div) <- paste(pheno_data_ENDO[1:24,1], pheno_data_ENDO[1:24,2], sep = "_" )

# color for shannon plot
col_shan <- c(unique(CoCuPp1), unique(CoCuSm1), unique(CoCuPp1), unique(CoCuSm1))

# Plot Shannon index
pdf(paste("ENDO_stats_plots/ms1_comp_list_diversity_shannon_Culture.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
#jpeg(filename = "ENDO_stats_plots/ms1_comp_list_diversity_shannon_Culture.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
boxplot(model_div$shannon ~ mzml_pheno_samples_type, col= col_shan, main="Shannon diversity (H\')", xlab="culture types", ylab="Shannon diversity index (H\')")
text(1:length(levels(mzml_pheno_samples_type)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples_type), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples_type))
text(1:length(levels(mzml_pheno_samples_type)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
legend("bottomleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(col_shan), legend= c("P.parvum", "S.marinoi"))
dev.off()

# ---------- PLS ----------
principal_components <- 5

# PLS according to culture type by species Pp
sel_pls_comp_list <- f.select_features_pls(feat_matrix=bina_list_Pp, 
                                           sel_factor=mzml_pheno_origin_samples[index_PP], 
                                           sel_colors=mzml_pheno_colors[index_PP], 
                                           components=principal_components, tune_length=10, 
                                           quantile_threshold=0.95, plot_roc_filename="ENDO_stats_plots/ms1_comp_list_select_pls_roc_Pp_bina.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=comp_list_Pp, 
                            sel_feat=sel_pls_comp_list$`_selected_variables_`, 
                            sel_names=paste0("",sel_pls_comp_list$`_selected_variables_`), 
                            sample_colors=mzml_pheno_colors[index_PP], plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename=NULL, main="PLS")
heatmaply(scale(comp_list_Pp[, which(colnames(comp_list_Pp) %in% sel_pls_comp_list$`_selected_variables_`)]), 
          k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
          file="ENDO_stats_plots/ms1_comp_list_select_pls_Pp_bina.html", 
          selfcontained=TRUE,
          column_text_angle = 90, 
          showticklabels = c(TRUE, TRUE),
          #titleX = FALSE,
          branches_lwd = 0.3,
          #width = 1000, heigth = 800,
          #colorbar_thickness = 50,
          #colorbar_len = 1,
          cexRow = 2, cexCol = 1.5, cexColorbar = 2)

sel_pls_comp_list$`_selected_variables_`
sel_pls_comp_list$`_model_r2_`
sel_pls_comp_list$`_multiclass_metrics_`

# selected features from PLS 
features_Pp <- sel_pls_comp_list$`_selected_variables_`
features_Pp <- sort(features_Pp, decreasing = FALSE)

#save(sel_pls_comp_list, file = "ENDO_stats_Results/sel_pls_comp_list_culture_Pp_bina.RData")
load("ENDO_stats_Results/sel_pls_comp_list_culture_Pp_bina.RData")
# which features in which condition P.PARVUM
bina_list_Pp_features <- bina_list_Pp[,features_Pp]
bina_list_Pp_features_Co <- bina_list_Pp_features[c(1,2,3,12),]
bina_list_Pp_features_Co <- bina_list_Pp_features_Co[,colSums(bina_list_Pp_features_Co)>0]
bina_features_Pp_Co <- data.frame(colnames(bina_list_Pp_features_Co))
bina_features_Pp_Co$Description <- "present in Pp Co-culture"
bina_features_Pp_Co$Condition <- "Co"
colnames(bina_features_Pp_Co) <- c("Feature", "Description","Condition")

bina_list_Pp_features_Mono <- bina_list_Pp_features[c(4:11),]
bina_list_Pp_features_Mono <- bina_list_Pp_features_Mono[,colSums(bina_list_Pp_features_Mono)>0]
bina_features_Pp_Mono <- data.frame(colnames(bina_list_Pp_features_Mono))
bina_features_Pp_Mono$Description <- "present in Pp Mono-culture"
bina_features_Pp_Mono$Condition <- "Mono"
colnames(bina_features_Pp_Mono) <- c("Feature", "Description","Condition")


sel_pls_features_Pp_bina <- data.frame(rbind(bina_features_Pp_Co, bina_features_Pp_Mono))
sel_pls_features_Pp_bina$Feature_id <- gsub("_neg","", sel_pls_features_Pp_bina$Feature)
sel_pls_features_Pp_bina$Feature_id <- gsub("_pos","", sel_pls_features_Pp_bina$Feature_id)
write.csv(sel_pls_features_Pp_bina, file = "ENDO_stats_Results/manuscript/sel_pls_features_Pp_bina.csv")



# PLS according to Culture type by species Sm
sel_pls_comp_list <- f.select_features_pls(feat_matrix=bina_list_Sm, 
                                           sel_factor=mzml_pheno_origin_samples[index_SM], 
                                           sel_colors=mzml_pheno_colors[index_SM], 
                                           components=principal_components, tune_length=10, 
                                           quantile_threshold=0.95, plot_roc_filename="ENDO_stats_plots/ms1_comp_list_select_pls_roc_Sm.pdf")
print(paste("Number of selected variables:", f.count.selected_features(sel_feat=sel_pls_comp_list$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=comp_list_Sm, 
                            sel_feat=sel_pls_comp_list$`_selected_variables_`, 
                            sel_names=paste0("",sel_pls_comp_list$`_selected_variables_`), 
                            sample_colors=mzml_pheno_colors[index_SM], plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename=NULL, main="PLS")
heatmaply(scale(comp_list_Sm[, which(colnames(comp_list_Sm) %in% sel_pls_comp_list$`_selected_variables_`)]), 
          k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
          file="ENDO_stats_plots/ms1_comp_list_select_pls_Sm_bina.html",
          column_text_angle = 90, 
          #width = 1000, heigth = 800, 
          showticklabels = c(TRUE, TRUE),
          #titleX = FALSE,
          branches_lwd = 0.3,
          
          selfcontained=TRUE, cexRow = 2, cexCol = 1.5, cexColorbar = 2)

sel_pls_comp_list$`_selected_variables_`
sel_pls_comp_list$`_model_r2_`
sel_pls_comp_list$`_multiclass_metrics_`

# selected features from PLS 
features_Sm <- sel_pls_comp_list$`_selected_variables_`
features_Sm <- sort(features_Sm, decreasing = FALSE)

save(sel_pls_comp_list, file = "ENDO_stats_Results/sel_pls_comp_list_culture_Sm_bina.RData")
load("ENDO_stats_Results/sel_pls_comp_list_culture_Sm_bina.RData")

features_Sm
# which features in which condition S.MARINOI
bina_list_Sm_features <- bina_list_Sm[,features_Sm]
bina_list_Sm_features_Co <- bina_list_Sm_features[c(1,2,3,12),]
bina_list_Sm_features_Co <- bina_list_Sm_features_Co[,colSums(bina_list_Sm_features_Co)>0]
bina_features_Sm_Co <- data.frame(colnames(bina_list_Sm_features_Co))
bina_features_Sm_Co$Description <- "present in Sm Co-culture"
bina_features_Sm_Co$Condition <- "Co"
colnames(bina_features_Sm_Co) <- c("Feature", "Description","Condition")

bina_list_Sm_features_Mono <- bina_list_Sm_features[c(4:11),]
bina_list_Sm_features_Mono <- bina_list_Sm_features_Mono[,colSums(bina_list_Sm_features_Mono)>0]
bina_features_Sm_Mono <- data.frame(colnames(bina_list_Sm_features_Mono))
bina_features_Sm_Mono$Description <- "present in Sm Mono-culture"
bina_features_Sm_Mono$Condition <- "Mono"
colnames(bina_features_Sm_Mono) <- c("Feature", "Description","Condition")


sel_pls_features_Sm_bina <- data.frame(rbind(bina_features_Sm_Co, bina_features_Sm_Mono))
sel_pls_features_Sm_bina$Feature_id <- gsub("_neg","", sel_pls_features_Sm_bina$Feature)
sel_pls_features_Sm_bina$Feature_id <- gsub("_pos","", sel_pls_features_Sm_bina$Feature_id)
write.csv(sel_pls_features_Sm_bina, file = "ENDO_stats_Results/manuscript/sel_pls_features_Sm_bina.csv")



# save selected features in a table
sel_pls_features_ENDO <- c(features_Pp, features_Sm)
origin_ENDO <- c(rep("features_Pp", length(features_Pp)), rep("features_Sm", length(features_Sm)))
feature <- gsub("_neg", "", sel_pls_features_ENDO)
feature <- gsub("_pos", "", feature)
sel_pls_feat_ENDO <- data.frame(origin_ENDO, sel_pls_features_ENDO, feature)
write.csv(sel_pls_feat_ENDO, file = "ENDO_stats_Results/sel_pls_feat_ENDO_bina.csv", row.names = FALSE)


# ---------- PLS on MS1 and MS2 object ----------
# load feature table MS1andMS2
feature_list_1ms2_pos <- read.csv("endo_pos/endo_pos_1ms2_Results/feature_list_endo_pos.csv", row.names = 1)
# remove MB
feature_list_1ms2_pos <- feature_list_1ms2_pos[-25,]

# set rownames of feature table
for (i in 1:24) {
  rownames(feature_list_1ms2_pos)[i] <- gsub("coculture_pos", pheno_data_ENDO$sample_group[i], rownames(feature_list_1ms2_pos)[i])
}
rownames(feature_list_1ms2_pos)[25:29] <- gsub("KSS_210401", "MS2", rownames(feature_list_1ms2_pos)[25:29])


comp_list_1ms2_pos <- feature_list_1ms2_pos[1:24,]

# where are features in MS2
test <- feature_list_1ms2_pos[25:29,]
test_feat <- which(colSums(test) != 0)
test_sub <- test[,test_feat]

# subset comp_list_1ms2 for features that contain MS2 info
comp_list_1ms2_pos_MS2 <- comp_list_1ms2_pos[,test_feat]


# ---------- Check CoCu Pp features in Sm samples ----------
# load bina list features
sel_pls_feat_ENDO_bina <- read.csv("ENDO_stats_Results/sel_pls_feat_ENDO_bina.csv")
features_Pp <- sel_pls_feat_ENDO_bina[1:30,2]

# index for monoculture and coculture for the species
index_MOSM <- grep("Mono-culture S.m", mzml_pheno$mzml_pheno_legend)
index_MOPP <- grep("Mono-culture P.p", mzml_pheno$mzml_pheno_legend)

# index for coculture for the species
index_COSM <- grep("Co-culture S.m", mzml_pheno$mzml_pheno_legend)
index_COPP <- grep("Co-culture P.p", mzml_pheno$mzml_pheno_legend)

# selected features from PLS for CoCuPp
length(features_Pp)

# search in Sm ENDO metabolome
comp_list_SmEn <- feat_list_ENDO[index_SM, which(colnames(comp_list_Sm) %in% features_Pp)]
comp_list_SmEn <- comp_list_SmEn[,colSums(comp_list_SmEn != 0) > 0]
# remove peaks that are only present in 1 sample
max(comp_list_SmEn)
comp_list_SmEn <- comp_list_SmEn[, (colSums(comp_list_SmEn) > max(comp_list_SmEn)) > 0]



# search in Monoculture Sm ENDO 
comp_list_MoSmEn <- comp_list[index_MOSM, which(colnames(comp_list_Sm) %in% features_Pp)]
comp_list_MoSmEn <- comp_list_MoSmEn[,colSums(comp_list_MoSmEn != 0) > 0]
# remove peaks that are only present in 1 sample
max(comp_list_MoSmEn)
comp_list_MoSmEn <- comp_list_MoSmEn[, (colSums(comp_list_MoSmEn) > max(comp_list_MoSmEn)) > 0]


# scale data
#comp_list_MoSmEn <- scale((comp_list[index_MOSM, which(colnames(comp_list_Sm) %in% features_CoCuPp)]))
#comp_list_MoSmEn[which(is.na(comp_list_MoSmEn))] <- 0

# heatmaply(comp_list_MoSmEn, 
#           k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
#           file="ENDO_stats_plots/ms1_comp_list_CoCuPp_in_Sm.html", selfcontained=TRUE)

# search in Coculture Sm ENDO
comp_list_CoSmEn <- comp_list[index_COSM, which(colnames(comp_list_Sm) %in% features_Pp)]
comp_list_CoSmEn <- comp_list_CoSmEn[,colSums(comp_list_CoSmEn != 0) > 0]
# remove peaks that are only present in 1 sample
max(comp_list_CoSmEn)
comp_list_CoSmEn <- comp_list_CoSmEn[, (colSums(comp_list_CoSmEn) > max(comp_list_CoSmEn)) > 0]

# heatmaply(comp_list_CoSmEn, 
#           k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
#           file="ENDO_stats_plots/ms1_comp_list_CoCuPp_in_SmCoCu.html", selfcontained=TRUE)
# 

# search in Monoculture Sm EXO 
load("exo_stats_Results/comp_list_EXO.RData")
comp_list_MoSmEx <- comp_list[index_MOSM, which(colnames(comp_list_Sm) %in% features_Pp)]
comp_list_MoSmEx <- comp_list_MoSmEx[,colSums(comp_list_MoSmEx != 0) > 0]
# remove peaks that are only present in 1 sample
max(comp_list_MoSmEx)
comp_list_MoSmEx <- comp_list_MoSmEx[, (colSums(comp_list_MoSmEx) > max(comp_list_MoSmEx)) > 0]
#subset again for which
comp_list_MoSmEx <- comp_list_MoSmEx[, which(colnames(comp_list_MoSmEx) %in% features_Pp)]



# plot in heatmap
# heatmaply(comp_list_MoSmEx, 
#           k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), 
#           file="ENDO_stats_plots/ms1_comp_list_CoCuPp_in_EXOSmMonoCu.html", selfcontained=TRUE)



# compare the found features
# ENDO Sm
features_ENDOSm_Pp <- colnames(comp_list_SmEn)
count_ENDOSm_Pp <- colSums(comp_list_SmEn != 0)
count_ENDOSm_Pp <- data.frame(count_ENDOSm_Pp)
nrow(count_ENDOSm_Pp)

write.csv(features_ENDOSm_Pp, "ENDO_stats_Results/manuscript/features_ENDOSm_Pp.csv")
# ENDO Sm Mono

# ENDO Sm Co


# EXO MonoCu Sm
features_EXOMonoSm_Pp <- colnames(comp_list_MoSmEx)
count_EXOMonoSm_Pp <- colSums(comp_list_MoSmEx != 0)
count_EXOMonoSm_Pp <- data.frame(count_EXOMonoSm_Pp)
nrow(count_EXOMonoSm_Pp)

write.csv(features_EXOMonoSm_Pp, "ENDO_stats_Results/manuscript/features_EXOMonoSm_Pp.csv")

# create one table with two rows that contain EXO Sm and ENDO Sm as column
ENDO_Sm <- NA
EXO_Sm <- NA
ENDO_Pp_feat <- data.frame(features_Pp, ENDO_Sm, EXO_Sm)
#rownames(ENDO_Pp_feat) <- rownames(count_EXOMonoSm_Pp)



library(VennDiagram)

plot(venn.diagram(
  x = list(Table1 = comp_list, Table2 = features_Pp),
  filename = "venn_diagram.png",  # Output file name
  col = "black",  # Set circle color
  fill = c("dodgerblue", "darkorange1"),  # Set fill color for each set
  alpha = 0.5,  # Set transparency level
  label.col = c("black", "black", "black"),  # Set label color for each set
  cex = 1.5,  # Set font size
  fontface = "bold",  # Set font face
  cat.col = c("dodgerblue", "darkorange1"),  # Set label color for intersection areas
  cat.cex = 1.5,  # Set font size for intersection labels
  cat.fontface = "bold",  # Set font face for intersection labels
  cat.pos = 0,  # Set label position for intersection areas
  cat.dist = 0.06,  # Set distance between label and area
  margin = 0.1  # Set margin size
))


# ---------- PCA  ----------
# combined feature table pos and neg
#pdf(paste("ENDO_stats_plots/PCA_legend.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=10, height=10, family="Helvetica")
jpeg(filename = "ENDO_stats_plots/ENDO_ms1_feature_table_pca.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
ms1_pca_ENDO <- prcomp(comp_list[1:24,], center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_ENDO$x[, 1], ms1_pca_ENDO$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color, cex=2)
grid()
text(ms1_pca_ENDO$x[,1], ms1_pca_ENDO$x[,2], labels=str_sub(rownames(feat_list_ENDO)[1:24], - 3, - 1), col=color, pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(mzml_pheno$mzml_pheno_legend))
dev.off()

# Skeletonema marinoi
jpeg(filename = "ENDO_stats_plots/ENDO_ms1_feature_table_pca_Sm.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
ms1_pca_ENDO <- prcomp(comp_list_Sm, center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_ENDO$x[, 1], ms1_pca_ENDO$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color[index_SM], cex=2)
grid()
text(ms1_pca_ENDO$x[,1], ms1_pca_ENDO$x[,2], labels=str_sub(rownames(feat_list_ENDO)[index_SM], - 3, - 1), col=color[index_SM], pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color[index_SM]), legend= unique(mzml_pheno$mzml_pheno_samples_type[index_SM]))
dev.off()

# Prymnesium parvum
jpeg(filename = "ENDO_stats_plots/ENDO_ms1_feature_table_pca_Pp.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
ms1_pca_ENDO <- prcomp(comp_list_Pp, center=TRUE)
par(mar=c(6,6,4,1), oma=c(0,0,0,0), cex.axis=2, cex=1, cex.lab=3, cex.main=3)
plot(ms1_pca_ENDO$x[, 1], ms1_pca_ENDO$x[,2], pch=19, main="PCA of feature table",
     xlab=paste0("PC1: ", format(summary(ms1_pca_ENDO)$importance[2, 1] * 100, digits=3), " % variance"),
     ylab=paste0("PC2: ", format(summary(ms1_pca_ENDO)$importance[2, 2] * 100, digits=3), " % variance"),
     col=color[index_PP], cex=2)
grid()
text(ms1_pca_ENDO$x[,1], ms1_pca_ENDO$x[,2], labels=str_sub(rownames(feat_list_ENDO)[index_PP], - 3, - 1), col=color[index_PP], pos=3, cex=1.5)
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color[index_PP]), legend= unique(mzml_pheno$mzml_pheno_samples_type[index_PP]))
dev.off()

# ---------- Variation partitioning ----------
# normalize comp_list
comp_list_varpart <- scale(comp_list_Pp)
# remove NaN from data 
comp_list_varpart[which(is.na(comp_list_varpart))] <- 0
# perfomr Varpart
model_varpart  <- varpart(comp_list_varpart, ~ mzml_pheno_origin_samples[index_PP] , ~ mzml_pheno_origin_samples[index_PP]  )

# Plot results
pdf(file="ENDO_stats_plots/pos_ms1_varpart_Pp.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart , Xnames=c("samples","culture type"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()


# --------- t-SNE ----------

library(Rtsne)

# ERROR: protection stack overflow when using the whole dataset
load("ENDO_stats_Results/comp_list_ENDO.RData")
tsne_ENDO <- Rtsne(test, perplexity = 50, check_duplicates = FALSE)
tsne_ENDO_Sm <- Rtsne(comp_list[index_SM,1:15000], perplexity = 2, check_duplicates = FALSE)
tsne_ENDO_Pp <- Rtsne(comp_list[index_PP,1:15000], perplexity = 2, check_duplicates = FALSE)

jpeg(filename = "ENDO_stats_plots/ENDO_ms1_tSNE_2.jpeg", width = 1000, height = 700, quality = 100, bg = "white")
par(mfrow=c(1,2), cex.axis=2, cex=1, cex.lab=1.5, cex.main=3)
plot(tsne_ENDO$Y, col = "black", bg= color, pch = 21, cex = 2,
     main="all conditions")
plot(tsne_ENDO_Sm$Y, col = "black", bg= color[index_SM], pch = 21, cex = 2,
     main="S. marinoi")
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color[index_SM]), legend= unique(mzml_pheno$mzml_pheno_samples_type[index_SM]))
#title(main = "t-SNE endometabolome")
plot(tsne_ENDO_Pp$Y, col = "black", bg= color[index_PP], pch = 21, cex = 2,
     main="P. parvum")
legend("topright", bty="n", pt.cex=5, cex=2, y.intersp=0.7,  inset = 0.2, xjust = 0, text.width=0.5, pch=20, 
       col= unique(color[index_PP]), legend= unique(mzml_pheno$mzml_pheno_samples_type[index_PP]))
#text(tsne_ENDO$Y[,1], tsne_ENDO$Y[,2], labels=str_sub(rownames(comp_list), - 3, - 1), col=color, pos=3, cex=1.5)
dev.off()


load("exo_stats_plots/comp_list_EXO.RData")
tsne_EXO <- Rtsne(comp_list[,1:15000], perplexity = 5, check_duplicates = FALSE)
tsne_EXO_Sm <- Rtsne(comp_list[index_SM,1:15000], perplexity = 2, check_duplicates = FALSE)
tsne_EXO_Pp <- Rtsne(comp_list[index_PP,1:15000], perplexity = 2, check_duplicates = FALSE)

jpeg(filename = "exo_stats_plots/EXO_ms1_tSNE.jpeg", width = 1500, height = 700, quality = 100, bg = "white")
par(mfrow=c(1,3), cex.axis=2, cex=1, cex.lab=1.5, cex.main=3)
plot(tsne_EXO$Y, col = "black", bg= color, pch = 21, cex = 2)
legend("topleft", bty="n", pt.cex=5, cex=2, y.intersp=0.7, text.width=0.5, pch=20, 
       col= unique(color), legend= unique(mzml_pheno$mzml_pheno_samples_type))
plot(tsne_EXO_Sm$Y, col = "black", bg= color[index_SM], pch = 21, cex = 2)
plot(tsne_EXO_Pp$Y, col = "black", bg= color[index_PP], pch = 21, cex = 2)
#text(tsne_EXO$Y[,1], tsne_EXO$Y[,2], labels=str_sub(rownames(comp_list), - 3, - 1), col=color, pos=3, cex=1.5)
dev.off()

