#  ---------------------SDM MODEL PREDICTIONS   -------------------------
# Species - African Tulip
# 1 load occurrence data
# 2 load climate/environmental data
# 3 prepare data for model
# 4 prepare biomod2 models (options and parameters set up)
# 5 run biomod2 model and project

# ------ Set Working Directory --------
setwd('C:/Data/sdm/C2_SDM/')
rm(list = ls()) # clean workspace

# install and load packages
# install.packages(c("dismo", "sdm", "fuzzySim", "maps", "terra", tidyverse", "ggplot2","sf", "enmSdmX"))
library(dismo)
library(fuzzySim)
library(maps)
library(ggplot2)
library(sf)
library(raster)
library(rJava)
options(java.parameters = "-Xmx4g")
library(terra)
library(tidyverse)
library(mecofun)

# ------- LOAD CLIM DATA ------
# load in raster datasets that have been saved (and have low correlation)
stck <- stack('./clim/cor_stack.tif')
pic_predictor <- stack('./clim/pic_stack.tif')

# update names of rasterstack
iv_names <- c("bio1", "bio2", "bio3", "bio7",
                           "bio12", "bio14", "bio18", "elevation")
names(stck) <- iv_names
names(pic_predictor) <- iv_names

# ------- LOAD OCCURRENCE DATA ------
# read in occurrence points that have been cleaned
occ = read.csv("C:/Data/sdm/C2_SDM/occurence/at/comb_at.csv", sep = ',')

# load in thinned data points (occurrence and absence data)
rec_sp_df = read.csv("C:/Data/sdm/C2_SDM/occurence/at/at_thinned_comb_all.csv", sep = ',')
pic_sp = read.csv("C:/Data/sdm/C2_SDM/occurence/at/at_thinned_comb_pic.csv", sep = ',')

# ----------------- PREPARE DATA FOR MODEL -----------------
# convert rec_sp_df to spatialpoints
rec_sp_xy <- SpatialPointsDataFrame(coords = rec_sp_df[,c(1,2)], data = rec_sp_df,
                                    proj4string = crs(stck))

# split all data into training and test dataset (mixed presence and absence)
rec_fold <- kfold(rec_sp_xy, k = 4)
rec_sp_test <- rec_sp_xy[rec_fold == 1,] # 1390
rec_sp_train <- rec_sp_xy[rec_fold != 1,] # 4168

# split thinned combined data based on presence and absence to generate
# individual abs and pres train/test datasets
pres_sp <- rec_sp_xy[rec_sp_xy$presence == 1,] # 2937
abs_sp <- rec_sp_xy[rec_sp_xy$presence == 0,] # 2621

# 25% sample for testing
fold <- kfold(pres_sp, k = 4)
occ_test <- pres_sp[fold == 1,] # 734
occ_train <- pres_sp[fold != 1,] # 2203

abs_fold = kfold(abs_sp, k = 4)
abs_test = abs_sp[abs_fold == 1,] # 655
abs_train = abs_sp[abs_fold != 1,] # 1966

# convert global train/test data to pacific
# get pacific region extent first
library(rnaturalearth)
world <- ne_download(type = "countries", scale = "medium")

# select only oceania region
pc <- subset(world, world$REGION_UN == "Oceania")
plot(pc)

# reproject pc first
pdc <- crs("+init=epsg:3832")
pc_rpj <- st_transform(pc, pdc)
plot(pc_rpj, max.plot = 1)

# get extent - ensure x object is spdf
pic_ext <- extent(pc_rpj)

pic_sp_train <- crop(spTransform(rec_sp_train, pdc), pic_ext) # 390
pic_sp_test <- crop(spTransform(rec_sp_test, pdc), pic_ext) # 148

# crop and reproject PRESENCE training/test data for PIC
pic_train <- crop(spTransform(occ_train, pdc), pic_ext) # 261
pic_test <- crop(spTransform(occ_test, pdc), pic_ext) # 92

# crop and reproject ABSENCE training/test data for PIC
pic_ab_train <- crop(spTransform(abs_train, pdc), pic_ext) # 144
pic_ab_test <- crop(spTransform(abs_test, pdc), pic_ext) # 41

# MESS calculation - global
test_df <- as.data.frame(rec_sp_test)
test_df <- test_df[test_df$presence == 1,]
test_df <- test_df[,c(4:11)] # 732 presence points extracted
dismo::mess(stck, test_df, filename = './results/mess_at_global.tif', overwrite = TRUE)

# MESS calculation - pacific region
pic_test_df <- as.data.frame(pic_test)
pic_test_df <- pic_test_df[pic_test_df$presence == 1,]
pic_test_df <- pic_test_df[,c(4:11)]
dismo::mess(pic_predictor, pic_test_df, filename = './results/mess_at_pic.tif', overwrite = TRUE)

# ------------------ Build Ensemble Model ------------------
library(sdm)
sdm_data <- sdmData(presence ~ .,
                    train = rec_sp_train[,c(3:11)],
                    test = rec_sp_test[,c(3:11)])
sdm_data

# fit using 3 models
model_cv <- sdm(presence~., data = sdm_data, methods = c('gam', 'rf', 'maxent'),
             replication = 'cv', cv.folds = 5, n = 1)

model_cv
# write sdm to file
write.sdm(model_cv, './results/sdm_model_cv_thinned', overwrite = T)
# model_cv <- read.sdm('C:/Data/sdm/C2_SDM/results/sdm_model_cv_thinned.sdm')

sdm_model_info = getModelInfo(model_cv)
write.csv(sdm_model_info, './results/sdm_m_info.csv')

# max(se+sp) opt for threshold-based stats (TSS)
sdm_model_eval = getEvaluation(model_cv, wtest = 'test.indep',
                               stat = c('AUC', 'COR','TSS', 'Deviance', 'threshold',
                                        'prevalence', 'sensitivity', 'specificity'),
                               opt = 2)
sdm_model_eval
write.csv(sdm_model_eval, './results/sdm_m_eval_all_thinned.csv')

roc(model_cv)

# get variable importance based on each model method
gam_vi = getVarImp(model_cv, method = 'gam')
plot(gam_vi, 'cor')

rf_vi = getVarImp(model_cv, method = 'rf')
plot(rf_vi, 'cor')

max_vi = getVarImp(model_cv, method = 'maxent')
plot(max_vi, 'cor')

gam_vi
rf_vi
max_vi

# get response curves for each model algorithm and run
rc = getResponseCurve(model_cv, id = 1:5) # response curve for GAM (model id 1-5)
rcurve(rc)
rc_df = rc@response
write.csv(rc_df, './results/rc_gam.csv')

rc_rf = getResponseCurve(model_cv, id = 6:10)
rcurve(rc_rf)
rc_rf_df = rc_rf@response
write.csv(rc_rf_df, './results/rc_rf.csv')

rc_max = getResponseCurve(model_cv, id = 11:15)
rcurve(rc_max)
rc_max_df = rc_max@response
write.csv(rc_max_df, './results/rc_maxent.csv')

# predict single models (to global and pacific (for comparison))
global_model_pred <- predict(model_cv, newdata = stck,
                      filename = "./results/global_pred_sdm_cv_thinned.tif", mean = T)

model_pred <- predict(model_cv, newdata = pic_predictor,
                      filename = "./results/predicted_sdm_cv_thinned.tif", overwrite = T)

# --------------- ENSEMBLE MODEL -------------------
# ensemble based on a Weighted averaging that is weighted using TSS statistic
en_model <- sdm::ensemble(model_cv, newdata = pic_predictor, filename = './results/ensemble_cv_thin2 .tif',
               setting = list(method = 'weighted',stat = 'TSS', opt = 2))
plot(en_model)

# get variable importance for ensemble
en_vi = getVarImp(model_cv, id = "ensemble", setting = list(method = 'weighted', stat = 'TSS', opt = 2))
en_vi

en_rc = getResponseCurve(model_cv, id = 'ensemble')
rcurve(en_rc)
en_rc_df = en_rc@response
write.csv(en_rc_df, './results/rc_ensembe.csv')

# evaluate model response based on training data
en_global_model <- sdm::ensemble(model_cv, newdata = as.data.frame(rec_sp_train[,c(3:11)]),
                                 filename = './results/en_global_cv_thin .tif',
                                 setting = list(method = 'weighted',stat = 'TSS', opt = 2))
en_global_model

en_e = evaluates(rec_sp_train$presence, en_global_model[,1])
en_e

# ------------- measure uncertainty --------------
# get uncertainty of model predictions
en_un_model <- sdm::ensemble(model_cv, newdata = pic_predictor, filename = './results/ensemble_uncertainty.tif',
               setting = list(method = 'uncertainty'))

# generate coefficient variation of probabilities generated from the multiple models
en_cv_model <- sdm::ensemble(model_cv, newdata = pic_predictor, filename = './results/ensemble_covar.tif',
               setting = list(method = 'cv'))

# --------------- Cont Boyce Index -------------------
# cont boyce index eval for pic prediction for GAM models
library(terra)
library(enmSdmX)

p_test = SpatialPoints(pic_test)
a_test = SpatialPoints(pic_ab_test)
model_all <- stack(model_pred)

# create function to calculate for pic region cbi for each model
cbi_sdm <- function(model){
  at_model_p <- terra::extract(model, p_test)
  at_model_a <- terra::extract(model, a_test)

  at_cbi <- evalContBoyce(pres = at_model_p,
                            contrast = at_model_a,
                            graph = T,
                            na.rm = T)
  at_cbi
}

# CBI score for GAM model
# gam
cbi_sdm(model_all$id_1__sp_presence__m_gam__re_cros) # 0.9233901
cbi_sdm(model_all$id_2__sp_presence__m_gam__re_cros) # 0.9059808
cbi_sdm(model_all$id_3__sp_presence__m_gam__re_cros) # 0.7135943
cbi_sdm(model_all$id_4__sp_presence__m_gam__re_cros) # 0.8829126
cbi_sdm(model_all$id_5__sp_presence__m_gam__re_cros) # 0.8691497

# RF
cbi_sdm(model_all$id_6__sp_presence__m_rf__re_cros) # 0.782165
cbi_sdm(model_all$id_7__sp_presence__m_rf__re_cros) # 0.8428235
cbi_sdm(model_all$id_8__sp_presence__m_rf__re_cros) # 0.8534443
cbi_sdm(model_all$id_9__sp_presence__m_rf__re_cros) # 0.9462
cbi_sdm(model_all$id_10__sp_presence__m_rf__re_cros) # 0.8454823

# Maxent

cbi_sdm(model_all$id_11__sp_presence__m_maxent__re_cros) # 0.9667797
cbi_sdm(model_all$id_12__sp_presence__m_maxent__re_cros) # 0.9527
cbi_sdm(model_all$id_13__sp_presence__m_maxent__re_cros) # 0.9729904
cbi_sdm(model_all$id_14__sp_presence__m_maxent__re_cros) # 0.9320251
cbi_sdm(model_all$id_15__sp_presence__m_maxent__re_cros) # 0.9283939

# ensemble
cbi_sdm(en_model) # 0.8662162