# ------ Set Working Directory --------
# All data saved in C drive
setwd('C:/Data/niche-modelling/SDM_22/devils/30sec')

rm(list=ls()) # clean workspace

# install and load packages
# install.packages(c("dismo","maptools","maps", "rJava", "tidyverse", 
#                    "ggplot2","raster","sp"))
library(dismo)
library(maptools)
library(maps)
library(ggplot2)
library(raster)
library(rJava)
library(sp)

library(tidyverse)
library(rgbif)
library(tidyverse)
library(CoordinateCleaner)

# ---------- DOWNLOAD SPECIES OCCURENCE DATA -------
# download species data
# GBIF Devils ivy - obtain updated name of scientific sp
devil_sp = occ_search(scientificName = "Epipremnum aureum", limit = 10000, hasCoordinate = TRUE) 

# select data from list
devil_data <- devil_sp$data

# add in ALA dataset - australian natural history museum
devil_csv <- read.csv('C:/Data/niche-modelling/SDM_22/devils/30sec/data/records-2023-03-13.csv')

# scale data down
devil_data <- devil_data %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)

# remove records without coordinates
devil_data <- devil_data %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

# identify duplicates based on lat and long coordinates
devil_dup = duplicated(devil_data[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 250 duplicates
length(devil_dup[devil_dup == TRUE])

# remove duplicates - 3711 remaining
devil_data <- devil_data[!devil_dup,]

# remove records without coordinates and duplicates for csv loaded
devil_csv <- devil_csv %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

devil_csv_dup <- duplicated(devil_csv[,c("decimalLongitude", "decimalLatitude")])
length(devil_csv_dup[devil_csv_dup == TRUE]) # get number of duplicates - 3 duplicates
devil_csv <- devil_csv[!devil_csv_dup,]

# scale data down - remove columns not useful (only select the following)
devil_csv <- devil_csv %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         dataResourceUid, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)


# convert to dataframe
devil_data <- data.frame(devil_data)

# plot using ggplot 
wrld = borders("world", colour = "gray50", fill = "gray50")

mplot <- ggplot() +  wrld +
  geom_point(data = devil_data, aes(x=decimalLongitude, y=decimalLatitude),
             colour = "red", size = 0.5) +
  geom_point(data = devil_csv, aes(x = decimalLongitude, y = decimalLatitude), colour = "blue", size = 1.5) +
  coord_sf()

mplot

# use coordinate cleaner to flag problem records
# install.packages('countrycode')
library(countrycode)
# convert country code ISO2 to ISO3
devil_data$countryCode <- countrycode(devil_data$countryCode, 
                                       origin = 'iso2c',
                                       destination = 'iso3c')

devil_csv$countryCode <- countrycode(devil_csv$countryCode, 
                                       origin = 'iso2c',
                                       destination = 'iso3c')

# flag coordinates not 'clean' based on tests
flags <- clean_coordinates(x = devil_data,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'outliers', 'seas', 'countries'))

# flag coordinates not 'clean' based on tests
flags_csv <- clean_coordinates(x = devil_csv,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'seas')) #countries

summary(flags) # 271 records flagged out of 1242
summary(flags_csv) # 73 flagged records
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
plot(flags_csv, lon = "decimalLongitude", lat = "decimalLatitude")

# exclude flagged records
devil_occ <- devil_data[flags$.summary,]
devil_occ_csv <- devil_csv[flags_csv$.summary,]

# merge datasets
devil_occ = subset(devil_occ, select=-c(gbifID))
devil_occ_csv = subset(devil_occ_csv, select = -c(dataResourceUid))

devil_data <- rbind(devil_occ, devil_occ_csv)

# identify duplicates based on lat and long coordinates
devil_occ_dup = duplicated(devil_data[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 3 duplicates (out of 827 observations)
length(devil_occ_dup[devil_occ_dup == TRUE])

# remove duplicates - returns 1991 observations
devil_data <- devil_data[!devil_occ_dup,]


# METADATA TEST
# remove records with low precision
hist(devil_data$coordinateUncertaintyInMeters / 1000, breaks = 20)

# removed coordinates with uncertainty up to 100km
devil_data <- devil_data %>%
  filter(coordinateUncertaintyInMeters/1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# identify basis of records - remove fossil, living & preserved specimen 
table(devil_data$basisOfRecord)

devil_data <- filter(devil_data, basisOfRecord == "HUMAN_OBSERVATION" |
                        basisOfRecord == "OBSERVATION" |
                        basisOfRecord == "OCCURENCE" |
                        basisOfRecord == 'PRESERVED_SPECIMEN') # kept due to low results from cleaning

table(devil_data$individualCount)

# remove individual counts that equate to 0 and higher than 100
devil_data <- devil_data %>%
  filter(individualCount > 0 | is.na(individualCount)) %>%
  filter(individualCount < 99 | is.na(individualCount))

# exclude old records- 1925 remove
table(devil_data$year)

# plot records
ggplot(devil_data, aes(x = year)) +
  geom_histogram(binwidth=1) +
  ggtitle("Year of Epiprenmum Occurences")

devil_data <- devil_data %>%
  filter(year > 1945)

# check family name
table(devil_data$family)

# check taxon rank
table(devil_data$taxonRank)

# separate data into Native Range
library(rgdal)
nr <- readOGR("C:/Data/niche-modelling/SDM_22/devils/devil_nr.shp")

# convert to df and add in columns (id and species)
nr_df <- as.data.frame(nr, row.names = T)
nr_df$id <- c("A")
nr_df$species <- c("Epipremnum aureum")

# convert native range to spdf (poly)
nr_shp <- SpatialPolygonsDataFrame(nr, nr_df, match.ID = TRUE)

# plot range
nrange <- fortify(nr_shp)

range <- ggplot() +
  wrld +
  geom_polygon(data = nrange, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title = element_blank())
range

# flag occurences within native range
range_flags <- cc_iucn(
  x = devil_data,
  range = nr_shp,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged"
)

# devil data within the range
data_fin <- devil_data[range_flags,]

range +
  geom_point(data = devil_data, aes(x = decimalLongitude, y = decimalLatitude), color = 'blue') +
  geom_point(data = data_fin, aes(x = decimalLongitude, y = decimalLatitude), color = 'red')

# ------------- Problematic datasets -------------
# ddmm to dd.dd conversion error - 0 records are flagged
out.ddmm <- cd_ddmm(devil_data, lon = "decimalLongitude",
                    lat = "decimalLatitude", 
                    ds = "species", diagnostic = T, 
                    diff = 1, value = "dataset")

# test for rasterized sampling
# cd_round function identifies datasets with significant proportion of coords
par(mfrow = c(2,2), mar = rep(2, 4))
out.round <- cd_round(data_fin, lon = "decimalLongitude", lat = "decimalLatitude",
                      ds = "species", value = "dataset",
                      T1 = 7, graphs = T)

# export occurence data as CSV for easier reference into working directory
library(sf)
st_write(devil_data, "devil_global.csv")
st_write(data_fin, "devil_nr.csv")

# ------- Access and Download Environmental Data ----------
# -------------- WORLDCLIM DATA -------------
# get WORLDCLIM data
library(geodata)
currentEnv = worldclim_global(var = "bio", res = 0.5, path = ".") # download 0.5 mins (30sec data)

# convert spatraster to raster
currentEnv.ras <- as(currentEnv, "Raster")
ext = extent(currentEnv.ras) #get extent

# rename rasterstack variables
names(currentEnv.ras) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",
                           "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", 
                           "bio16", "bio17", "bio18", "bio19")
# ------------ DROP BIOCLIMATIC VARIABLES -----------------------------
# Limit set of bioclimatic predictors
# available list of bioclimatic indicators - https://www.worldclim.org/data/bioclim.html 
currentEnv.ras = dropLayer(currentEnv.ras, c("bio2", "bio3", "bio4", "bio10","bio11", 
                                     "bio13", "bio14", "bio15", "bio18", "bio19"))

# write each raster to file
writeRaster(currentEnv.ras, filename = names(currentEnv.ras), bylayer = TRUE, format = 'GTiff', overwrite=TRUE)

# ------- ELEVATION VARIABLES ---------- #
# obtain other environmental variables - elevation, slope, hillshade
# download elevation variables
elevation_wrld <- worldclim_global(var="elev", res=0.5, path = ".") # alt from srtm
#convert to raster
elevation_wrld <- as(elevation_wrld, "Raster")
raster::writeRaster(elevation_wrld, filename='world_elevation.tif', format='GTiff', overwrite=TRUE)

# plot paragraphs
par(mfrow = c(1,1))
plot(elevation_wrld)

# run topographic analysis
# calculate slope
slp = terra::terrain(elevation_wrld, v = 'slope', unit = 'radians', neighbors = 8, 
                     filename = 'slope_wrld.tif', overwrite = TRUE)
# convert spat raster to raster layer
slp <- as(slp, "Raster")
plot(slp)

# -------- Categorical Variable - Terrestrial Ecoregions --------
# read in raster biome
biome <- raster("C:/Data/niche-modelling/SDM_22/merremia/data/biome/wrld_biome.tif")

# clip biome to extent - remove polar regions (south)
ext_biome <- crop(biome, ext)

# resample biome layer to match extent, crs and resolution of global elevation before stack
stck_biome <- resample(ext_biome, elevation_wrld, method="bilinear")
writeRaster(stck_biome, filename="biome.tif", format='GTiff', overwrite=TRUE)

# --------MERGE DOWNLOADED VARIABLES----------
# stack raster layers
datafiles = Sys.glob("*.tif")
stck = stack()
for (i in 1:NROW(datafiles)){
  tempraster = raster(datafiles[i]) #iterate each file and store as temporary raster
  stck = stack(stck, tempraster)
}

# get list of names from file
r.ls <- (list.files(pattern = "tif$"))
# rename stack with layer names
names(stck) <- r.ls
# check names
names(stck)

# --------------------Check CORRELATION ------------------
library(sdmpredictors)

# run correlation co-efficient
p_mat <- pearson_correlation_matrix(stck, same_mask = FALSE)
View(p_mat) # remove layers with high correlation > 0.8 

# drop highly coorelated predicts
stck = dropLayer(stck, c("bio1.tif", "bio6.tif", "bio9.tif", "bio17.tif"))

# ------- GENERATE NATIVE RANGE PREDICTOR --------------------
# generate extent from native range
ext_nr <- raster::extent(nr_shp)
nr_predictor <- crop(stck, ext_nr, snap = 'near', datatype = NULL) # outputs raster brick
names(nr_predictor)

plot(nr_predictor)

# --------------------- PARTITION DATA ------------------
#------------------ SPLIT PRESENCE POINTS----------------
# 20% sample for testing, 80% for training
# global data
devil_data_wrld <- devil_data[,2:3]
devil_data_wrld <- SpatialPoints(devil_data_wrld, crs(stck))
fold_wrld <- kfold(devil_data_wrld, k = 5)
devil_test_wrld <- devil_data_wrld[fold_wrld == 1,]
devil_train_wrld <- devil_data_wrld[fold_wrld != 1,]

# native range data
nr_data <- data_fin[,2:3]
nr_data <- SpatialPoints(nr_data, crs(stck))
fold <- kfold(nr_data, k = 5)
nr_test <- nr_data[fold == 1,] # 20% of total devil data from native range
nr_train <- nr_data[fold != 1,] # 80 of total devil data from native range

# ---------------- ABSENCE POINTS --------------
# create background points/pseudo-absences for global and native range
bck = dismo::randomPoints(stck, n = 1500, p = devil_data_wrld, ext = extent(stck), extf = 1.1)
nr_bck = dismo::randomPoints(nr_predictor, n = 500, p = nr_data, ext = ext_nr, extf = 1.1)

colnames(bck) = c('lon', 'lat') # name columns for background points
colnames(nr_bck) = c('lon', 'lat')

# create sp layer with coordinate reference system
bck <- SpatialPoints(bck, crs(stck))
nr_bck <- SpatialPoints(nr_bck, crs(nr_predictor))

plot(bck, pch = 20) # check background points
plot(nr_bck, pch = 20, col = 'red', add = TRUE)

# Partition Data 
# group background points/pseduo-absences - for global range
group = kfold(bck, k = 5)
# group pseudo-absences for training model performances
bck_test = bck[group == 1,]
bck_train = bck[group != 1,]

# group background points/pseduo-absences - for native range
nr_group = kfold(nr_bck, k = 5)
# group pseudo-absences for training model performances
nr_bck_test = nr_bck[nr_group == 1,]
nr_bck_train = nr_bck[nr_group != 1,]

#------------- BUILD MAXENT ----------------#
# fit maxent - using stacked raster predictor variables
# PRESENCE ONLY
# global model with global distribution
maxent_model_wrld <- maxent(stck, devil_train_wrld, factor = "biome.tif") 
maxent_model_wrld

# native range model (predictor and test)
maxent_model_nr <- maxent(nr_predictor, nr_train, factor = "biome.tif")
maxent_model_nr

# plot and view variable contribution
par(mfrow = c(1,2))
plot(maxent_model_wrld) 
plot(maxent_model_nr) 

# run response curve for each model and variable
response(maxent_model_wrld)

response(maxent_model_nr)

# drop predictors with that contribute no response
# stck = dropLayer(stck, c("bio8.tif", "bio5.tif"))
# nr_predictor = dropLayer(nr_predictor, c("bio8.tif", "bio5.tif"))

# -------------- EVALUATE RUNNING AUC ------------------
# global prediction (on presence only)
e_wrld = evaluate(devil_test_wrld, bck_test, maxent_model_wrld, stck)

# e = evaluate(nr_test, nr_bck_test, maxent_model, stck)
e_nr =evaluate(nr_test, nr_bck_test, maxent_model_nr, nr_predictor)

e_wrld
e_nr

# plot roc and tpr curve
par(mfrow = c(1,2))
plot(e_wrld, 'ROC', main = 'ROC Curve - Maxent Current Predictor Variables')
plot(e_wrld, 'TPR', main = 'TPR Curve - Maxent Current Preditor Variables')

par(mfrow = c(1,2))
plot(e_nr, 'ROC', main = 'ROC Curve - Maxent Current Predictor Variables')
plot(e_nr, 'TPR', main = 'TPR Curve - Maxent Current Preditor Variables')

# ------------------ RUN PREDICTION ON LARGER AREA ---------------- 
# --------------------- CROP DATA TO PACIFIC REGION ------------
# restrict data to pacific/oceania region (def as region 5)
data(wrld_simpl)
pc <- subset(wrld_simpl, wrld_simpl$REGION==9)
# reproject pc first
pdc <- crs.latlong <- crs("+init=epsg:3832")
pc_rpj <- spTransform(pc, crs(pdc))
plot(pc_rpj)

# get extent
pic_ext <- extent(pc_rpj)

# ------- CROP PREDICTOR VARIABLES --------
# get predictor variables for PIC region
pic_predictor <- projectRaster(stck, crs = pdc) # fairly long processing - alternative to reproject via GIS software

# crop reprojected rasters based on extent of pacific region shp
pic_predictor <- crop(pic_predictor, pic_ext)

# ---------------------------- RUN PREDICTION ON PACIFIC REGION -----------------------------
# Presence Only Models
# predict the entire training dataset and export prediction results as raster
devil_pred_wrld <- predict(maxent_model_wrld, pic_predictor, filename = 'devil_predicted_global.tif', ext = pic_ext)

devil_pred_nr <- predict(maxent_model_nr, pic_predictor, filename = 'devil_predicted_nr.tif',ext = pic_ext)

#------- COLOR RAMP --------
# setting color ramp
clr <- RColorBrewer::brewer.pal(9, "YlOrRd")
rdcl <- RColorBrewer::brewer.pal(9, 'Reds')
crmp <- colorRampPalette(clr)
rdcrmp <- colorRampPalette(rdcl)
par(mfrow=c(1,1))

# plot predicted models
# pres only models
par(mfrow=c(1,2))
plot(devil_pred_wrld, col = rdcrmp(10), main = 'Devils Ivy Predicted Suitability under Current Global Model')
maps::map(database ='world', fill=FALSE, add = TRUE)

plot(devil_pred_nr, col = rdcrmp(10), main ='Devils Ivy Predicted Suitability under Current Native Range Model')
maps::map(database = 'world', fill = FALSE, add = TRUE)


# --------------------- BOYCE INDEX --------------------#
# run boyce index for model evaluation
remotes::install_github('adamlilith/enmSdm', dependencies=TRUE)
library(enmSdm)
library(terra)

# split train and test occurrence datasets again
devil_data_pic <- spTransform(devil_data_wrld, crs(pdc))
devil_data_pic <- crop(devil_data_pic, pic_ext)

# split training and test (80/20)
fold_pic <- kfold(devil_data_pic, k = 5)
devil_pic_test <- devil_data_pic[fold_pic == 1,]
devil_pic_train <- devil_data_pic[fold_pic != 1,]

# generate PIC background absence points 
devil_pic_bck_test <- spTransform(bck_test, crs(pdc))
devil_pic_bck_test <- crop(devil_pic_bck_test, pic_ext)

devil_pic_bck_train <- spTransform(bck_train, crs(pdc))
devil_pic_bck_train <- crop(devil_pic_bck_train, pic_ext)


# --------------PRESENCE ONLY MODELS
# global predictors and presence points
nr_prob_1 <- terra::extract(devil_pred_wrld, devil_pic_test)
prob_abs_1 <- terra::extract(devil_pred_wrld, devil_pic_bck_test)

contBoyce(pres = nr_prob_1, 
          contrast = prob_abs_1,
          na.rm = T)

# native predictors and native range presence points
nr_prob_1 <- terra::extract(devil_pred_nr, devil_pic_test)
prob_abs_1 <- terra::extract(devil_pred_nr, devil_pic_bck_test)

contBoyce(pres = nr_prob_1, 
          contrast = prob_abs_1,
          na.rm = T)
