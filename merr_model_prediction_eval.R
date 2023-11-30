# ------ Set Working Directory --------
setwd('C:/Data/niche-modelling/SDM_22/merremia/30sec/data')
rm(list=ls()) # clean workspace

# install and load packages
install.packages(c("dismo", "rgbif", "maps", "rJava", "tidyverse", "ggplot2", "raster","sp", "CoordinateCleaner"))
library(dismo)
library(maps)
library(ggplot2)

library(sp)
library(raster)
library(rJava)

library(rgbif)
library(terra)
library(tidyverse)
library(CoordinateCleaner)

# read in occurrnence points that have been cleaned
merr_occ = read.csv("merr_global_1.csv", sep = ',')
data_fin = read.csv("merr_nr_1.csv", sep = ',')

# ----- PARTITION DATA ---------
# PRESENCE POINTS
# 20% sample for testing
merr_df = cbind.data.frame(merr_occ$decimalLongitude, merr_occ$decimalLatitude)
merr_df <- SpatialPoints(merr_df, crs(stck))
fold <- kfold(merr_df, k=5)
merr_test <- merr_df[fold == 1,]
merr_train <- merr_df[fold!=1,]

# native range data 
nr_data <- data_fin[,2:3]
nr_data <- SpatialPoints(nr_data, crs(stck))
fold <- kfold(nr_data, k = 5)
nr_test <- nr_data[fold == 1,] # 20% of total occ data from native range
nr_train <- nr_data[fold != 1,] # 80 of total occ data from native range

# ABSENCE POINTS
# generate background points/pseudo-abs
ext = extent(stck) #get extent
merr_bck = randomPoints(stck, n = 800,p = merr_df, ext = ext, extf = 1.1)
nr_bck = randomPoints(nr_predictor, n = 700, p = nr_data, ext = ext_nr, extf = 1.1)

# column names for both background_points
colnames(merr_bck) = c('lon', 'lat')
colnames(nr_bck) = c('lon', 'lat')

# create sp layer with coordinate reference system for global generated bckgrd points
merr_bck <- SpatialPoints(merr_bck, crs(slp))
merr_group = kfold(merr_bck, k = 5)
merr_bck_test = merr_bck[merr_group == 1,]
merr_bck_train = merr_bck[merr_group != 1,]

merr_nr_bck <- SpatialPoints(nr_bck, crs(nr_predictor))
merr_nr_bck_group = kfold(merr_nr_bck, k = )
merr_nr_bck_test = merr_nr_bck[merr_nr_bck_group == 1,]
merr_nr_bck_train = merr_nr_bck[merr_nr_bck_group != 1,]

plot(merr_bck, pch = 20, col = 'blue') # check globally generated background points
plot(merr_nr_bck, pch = 20, col = 'red', add = TRUE)


# ----------------LOAD IN PREDICTORS ---------------------------------
# ------------MERGE DOWNLOADED VARIABLES---------------
# load in raster datasets that have been saved
datafiles = Sys.glob("*.tif")
stck = stack()

for (i in 1:NROW(datafiles)){
  tempraster = raster(datafiles[i])#iterate each file and store as temporary raster
  stck = stack(stck, tempraster)
}

# get list of names from file
r.ls <- (list.files(pattern = "tif$"))
# rename stack with layer names
names(stck) <- r.ls
# check names
names(stck)

# ------------------RUN CORRELATION------------------
install.packages('sdmpredictors')
library(sdmpredictors)

# run correlation co-efficient
p_mat <- pearson_correlation_matrix(stck, same_mask = FALSE)
View(p_mat) # remove layers with high correlation > 0.8 (bio 12 and bio 16)

# drop highly correlated predicts
stck = dropLayer(stck, c("bio1.tif", "bio6.tif", "bio17.tif", "world_aspect.tif", 
                         "preciptn.tif", "tempavg.tif", "tempmin.tif", "bio13.tif"))


# ------- GENERATE NATIVE RANGE PREDICTOR --------------------
# generate extent from native range
# Split occurence data within defined Native Range
nr <- vect("C:/Data/niche-modelling/SDM_22/merremia/30sec/nr_merr.shp")
nr <- as(nr, "Spatial")
# convert to df and add in columns (id and species)
nr_df <- as.data.frame(nr, row.names = T)
nr_df$id <- c("A")
nr_df$species <- c("Decalobanthus peltatus")

# convert native range extent to a spdf (poly)
nr_shp <- SpatialPolygonsDataFrame(nr, nr_df, match.ID = TRUE)
nr <- vect(nr_shp)
ext_nr <- raster::extent(nr_shp)
nr_predictor <- crop(stck, ext_nr, snap = 'near', datatype = NULL) # outputs raster brick
names(nr_predictor)

plot(nr_predictor)

#---------- FIT MAXENT ------------#
# fit maxent - using stacked raster predictor variables
# PRES ONLY MODELS
 # Global range model
model <- maxent(stck, merr_train, factor = "biome.tif")
# NRP
# model_nrp <- maxent(stck, nr_train, factor = "biome.tif")
# Native range model
nr_model <- maxent(nr_predictor, nr_train, factor = "biome.tif")

par(mfrow = c(1,2))
plot(model, main = "Presence Only Global Merremia Maxent Model")
#plot(model_nrp, main = "Pres Only Native Range Occurrence Model")
plot(nr_model, main = "Pres Only Native Range Restricted Model")

# run response curve
# global response model
response(model, expand = 0, rug = FALSE)
# native range response model
response(nr_model, expand = 0, rug = FALSE)

# remove low performing variables and then re-run model model!!!!!
stck = dropLayer(stck, c("bio5.tif", "bio12.tif", "bio7.tif"))
nr_predictor = dropLayer(nr_predictor, c("bio5.tif", "bio12.tif", "bio7.tif"))

# ---------- EVALUATE TRAINED MODEL -------------#
# run evaluation - using predictor variables
e1 = evaluate(merr_test, merr_bck_test, model, stck) 
e1

e2 = evaluate(nr_test, merr_nr_bck_test, nr_model, nr_predictor) 
e2

# plot ROC - and AUC
par(mfrow = c(1,2))
plot(e1, 'ROC', main = 'ROC Curve - Maxent Current Predictor Variables')
plot(e1, 'TPR', main = 'TPR Curve - Maxent Current Preditor Variables')

plot(e2, 'ROC', main = 'ROC Curve - Maxent Current Predictor Variables')
plot(e2, 'TPR', main = 'TPR Curve - Maxent Current Preditor Variables')

# ------------------ RUN PREDICTION ON LARGER AREA ---------------- 
# --------------------- CROP DATA TO PACIFIC REGION ------------
# restrict data to pacific/Oceania region 
library(rnaturalearth)
world <- ne_download(type = "countries", scale = "medium")

# select only oceania region
pc <- subset(world, oceania$REGION_UN=="Oceania")
plot(pc)

# re project pc first
pdc <- crs.latlong <- crs("+init=epsg:3832")
pc_rpj <- spTransform(pc, crs(pdc))
plot(pc_rpj)

# get extent
pic_ext <- extent(pc_rpj)

# ------------------- CROP PREDICTOR VARIABLES ------------------
# get predictor variables for PIC region - already reprojected to PDC
setwd('C:/Data/niche-modelling/SDM_22/merremia/30sec/pic_crop2')
pic_datafiles = Sys.glob("*.tif")
pic_predictor = stack()

for (i in 1:NROW(pic_datafiles)){
  tempraster = raster(pic_datafiles[i]) #iterate each file and store as temporary raster
  tempraster <- crop(tempraster, pic_ext)
  pic_predictor = stack(pic_predictor, tempraster)
}

# get list of names from file
r.ls <- (list.files(pattern = "tif$"))

# rename stack with layer names
names(pic_predictor) <- r.ls
# check names
names(pic_predictor)

# --------------------- RUN PREDICTION -------------------------------
# predict  trained models - pres only
merr_predict <- predict(model, pic_predictor, ext = pic_ext)
merr_predict_nr <- predict(nr_model, pic_predictor, ext = pic_ext)

#--------------------- COLOR RAMP ------------------
# setting color ramp
clr <- RColorBrewer::brewer.pal(9, "YlOrRd")
rdcl <- RColorBrewer::brewer.pal(9, 'Reds')
crmp <- colorRampPalette(clr)
rdcrmp <- colorRampPalette(rdcl)

#par(mfrow=c(2,3))
par(mfrow=c(1,2))
# plot predicted model - presence only models
plot(merr_predict, col = rdcrmp(10), main = 'Global Predicted Suitability for Merremia - Pres Only')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_predict, 'merremia_predicted_global2.tif')

plot(merr_predict_nr, col = rdcrmp(10), main = 'Native Range Predicted Suitability for Merremia - Pres Only')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_predict_nr, 'merremia_predicted_nr2.tif')


# EVALUATION
# ----------------------------- BOYCE INDEX --------------------------#
# run boyce index for model evaluation
# remotes::install_github('adamlilith/enmSdm', dependencies=TRUE)
library(enmSdm)

# CREATE PRESENCE POINT DATA
merr_df_pic <- spTransform(merr_df, crs(pdc))
merr_df_pic <- crop(merr_df_pic, pic_ext) 

# CREATE ABSENCE POINT DATA ---------  
# create random points within the pacific region (absence points) 
# generated random points within the pacific region extent
merr_bck_pic = randomPoints(pic_predictor, n = 150, p = merr_df_pic, ext = pic_ext, extf = 1.1)
merr_bck_pic = SpatialPoints(merr_bck_pic, crs(pic_predictor))

# -------- BOYCE INDEX
# PRESENCE  MODEL PREDICTION
# Pres - Merremia Probability for Pacific Region
# Global presence only
merr_prob_pres <- raster::extract(merr_predict, merr_df_pic)
# crop testing background absence points
merr_prob_abs <- raster::extract(merr_predict, merr_bck_pic)
merr_gb_cont <- contBoyce(pres = merr_prob_pres, 
          contrast = merr_prob_abs,
          na.rm = T)
merr_gb_cont

# native range presence and predictor
merr_prob_pres_nr <- raster::extract(merr_predict_nr, merr_df_pic)
merr_prob_abs_nr <- raster::extract(merr_predict_nr, merr_bck_pic)
merr_nr_cont <- contBoyce(pres = merr_prob_pres_nr, 
          contrast = merr_prob_abs_nr,
          na.rm = T)
merr_nr_cont