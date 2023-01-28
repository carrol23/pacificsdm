# reference material 1. - https://cmerow.github.io/RDataScience/3_6_Teaching_Ecoinformatics.html
# reference material 2 - https://consbiol-unibern.github.io/SDMtune/articles/articles/prepare-data.html
# reference material 3 - coordinate cleaner : https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_GBIF_data_with_CoordinateCleaner.html
# global human footprint - https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint 

# INSTALL NECESSARY PACKAGES
#install.packages("devtools")
library(devtools)
## install_github("ropensci/CoordinateCleaner")

library(CoordinateCleaner)

# load required libraries
library(rgbif)
# library(rgdal)
library(maps)
library(dismo) # ensure loaded well - maxent modelling package
library(maptools)
library(ggplot2)

# install.packages('geodata')
library(geodata)
library(rJava)

library(jsonlite)

library(caret) # loads ggplot2 and lattice as required

# ------ Set Working Directory --------
setwd('C:/Data/niche-modelling/SDM_22/spathodea/30sec_model')

rm(list=ls()) # clean workspace

# ------------ Download Species Data -------------
# obtain using rgbif
atulip_data = occ_search(scientificName = "Spathodea campanulata", limit = 10000, hasCoordinate = TRUE)
# select data from list
atulip_data <- atulip_data$data

library(dplyr)
# scale data down
atulip_data <- atulip_data %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)

# remove records without coordinates
atulip_data <- atulip_data %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# identify duplicates based on lat and long coordinates
atulip_d = duplicated(atulip_data[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 802 duplicates
length(atulip_d[atulip_d == TRUE])

# remove duplicates - 3711 remaining
atulip_data <- atulip_data[!atulip_d,]

# convert to dataframe
atulip_data <- data.frame(atulip_data)

# plot using ggplot 
wrld = borders("world", colour = "gray50", fill = "gray50")

ggplot() +  wrld +
  geom_point(data = atulip_data, aes(x=decimalLongitude, y=decimalLatitude),
             colour = "red", size = 0.5) +
  coord_sf()

# use coordinate cleaner to flag problem records
# install.packages('countrycode')
library(countrycode)
# convert country code ISO2 to ISO3
atulip_data$countryCode <- countrycode(atulip_data$countryCode, 
                                       origin = 'iso2c',
                                       destination = 'iso3c')

# flag coordinates not 'clean' based on tests
flags <- clean_coordinates(x = atulip_data,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'outliers', 'seas', 'countries'))

summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# exclude flagged records
atulip_data <- atulip_data[flags$.summary,]
# flagged records
atulip_flag <- atulip_data[!flags$.summary,]

# use metadata to remove data
# remove records with low precision
hist(atulip_data$coordinateUncertaintyInMeters / 1000, breaks = 20)

# removed coordinates with uncertainty up to 100km
atulip_data <- atulip_data %>%
  filter(coordinateUncertaintyInMeters/1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# identify basis of records - remove fossil, living & preserved specimen 
table(atulip_data$basisOfRecord)

atulip_data <- filter(atulip_data, basisOfRecord == "HUMAN_OBSERVATION" |
                        basisOfRecord == "OBSERVATION" |
                        basisOfRecord == "OCCURENCE")

table(atulip_data$individualCount)

# remove individual counts that equate to 0 and higher than 100
atulip_data <- atulip_data %>%
  filter(individualCount > 0 | is.na(individualCount)) %>%
  filter(individualCount < 99 | is.na(individualCount))

# exclude old records- 1925 remove
table(atulip_data$year)

# plot records
ggplot(atulip_data, aes(x = year)) +
  geom_histogram(binwidth=1) +
  ggtitle("Year of African Tulip Occurences")

atulip_data <- atulip_data %>%
  filter(year > 1945)

# check family name
table(atulip_data$family)

# check taxon rank
table(atulip_data$taxonRank)

# ------------------- NATIVE RANGE EXTENT -----------------
# import generated extent of atulip native range
library(rgdal)
nr <- readOGR("./shp", "nr_ext")
# convert to df and add in columns (id and species)
nr_df <- as.data.frame(nr, row.names = T)
nr_df$id <- c("A")
nr_df$species <- c("Spathodea campanulata")

# convert native range to spdf (poly)
nr_shp <- SpatialPolygonsDataFrame(nr, nr_df, match.ID = TRUE)

# plot range
nrange <- fortify(nr_shp)

at_range <- ggplot() +
  wrld +
  geom_polygon(data = nrange, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title = element_blank())

# show plot
at_range

# flag occurences within native range
range_flags <- cc_iucn(
  x = atulip_data,
  range = nr_shp,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged"
)

# atulip data within the range
data_fin <- atulip_data[range_flags,]

at_range +
  geom_point(data = atulip_data, aes(x=decimalLongitude, y = decimalLatitude), color = 'blue') +
  geom_point(data = data_fin, aes(x=decimalLongitude, y=decimalLatitude), color = 'red')

# ------------- Problematic datasets -------------
# ddmm to dd.dd conversion error - 0 records
out.ddmm <- cd_ddmm(atulip_data, lon = "decimalLongitude",
                    lat = "decimalLatitude", 
                    ds = "species", diagnostic = T, 
                    diff = 1, value = "dataset")

# test for rasterized sampling
# cd_round function identifies datasets with significant proportion of coords
par(mfrow = c(2,2), mar = rep(2, 4))
out.round <- cd_round(data_fin, lon = "decimalLongitude", lat = "decimalLatitude",
                      ds = "species", value = "dataset",
                      T1 = 7, graphs = T)

# ------- Access and Download Environmental Data ----------
# -------------- WORLDCLIM DATA -------------
# get WORLDCLIM data
library(geodata)
currentEnv = worldclim_global(var = "bio", res = 0.5, path = ".") # download 0.5 mins (30sec data0)

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
currentEnv.ras = dropLayer(currentEnv.ras, c("bio1", "bio3", "bio4", "bio8",
                                     "bio9", "bio10", "bio15", 
                                     "bio16", "bio17", "bio18", "bio19"))

# write each raster to file
writeRaster(stack(currentEnv.ras), names(currentEnv.ras), bylayer = TRUE, format = 'GTiff', overwrite=TRUE)

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
biome <- raster("./biome/wrld_biome.tif")
# clip biome to extent - remove polar regions (south)
ext_biome <- crop(biome, ext)

# resample biome layer to match extent, crs and resolution of elevation_wrld
stck_biome <- resample(ext_biome, elevation_wrld, method="ngb")

# assign projection from slp layer to biome dataset
stck_biome <- projectRaster(stck_biome, crs = crs(slp))
writeRaster(stck_biome, filename="biome.tif", format='GTiff', overwrite=TRUE)

# ---------------MERGE DOWNLOADED PREDICTOR VARIABLES-----------------
# stack raster layers
datafiles = Sys.glob("*.tif")
stck = stack()

for (i in 1:NROW(datafiles)){
  tempraster = raster(datafiles[i]) #iterate each file and store as temporary raster
  tempraster = (datafiles[[i]])
  stck = stack(stck, tempraster)
}

# get list of names from file
r.ls <- (list.files(pattern = "tif$"))
# rename stack with layer names
names(stck) <- r.ls
# check names
names(stck)
# -------------- run correlation co-efficient -----------------
# install sdmpredictor package
# install.packages('sdmpredictors')
library(sdmpredictors)

# run correlation co-efficient
p_mat <- pearson_correlation_matrix(stck, same_mask = FALSE)
View(p_mat) # remove layers with high correlation > 0.8 (bio 12 and bio 16)

# ensure all raster layers are cropped to the same extent
predictor <- crop(stck, ext)

# drop highly coorelated predicts
predictor = dropLayer(predictor, c("bio11.tif", "bio13.tif"))

# generate extent from native range
ext_nr <- extent(nr_shp)
nr_predictor <- crop(predictor, ext_nr, snap = 'near', datatype = NULL) # outputs raster brick
names(nr_predictor)

# check predictors
plot(predictor)

# --------------------- PARTITION DATA ------------------
#------------------ SPLIT PRESENCE POINTS----------------
# 20% sample for testing
# global data

atulip_data_wrld <- atulip_data[,2:3]
atulip_data_wrld <- SpatialPoints(atulip_data_wrld, crs(slp))
fold_wrld <- kfold(atulip_data_wrld, k = 5)
atulip_test_wrld <- atulip_data_wrld[fold_wrld == 1,]
atulip_train_wrld <- atulip_data_wrld[fold_wrld!=1,]

# native range data
nr_data <- data_fin[,2:3]
nr_data <- SpatialPoints(nr_data, crs(slp))
fold <- kfold(nr_data, k=5)
atulip_test <- nr_data[fold == 1,] # 20% of total atulip data from native range
atulip_train <- nr_data[fold!=1,] # 80 of total atulip data from native range

# ------------------ABSENCE POINTS-------------------
# create background points/pseudo-absences - p must be in a spatialdatapoint
bck = randomPoints(predictor, n=3000, p=atulip_data_wrld , ext=ext, extf = 1.1)

# crop existing raster extent to the native range 
nr_slp <- crop(slp, ext_nr)
# generate background points/pseudo-abs within the native range (use presence points constrainted within native range)
nr_bck = randomPoints(nr_slp, n=300, p=nr_data, ext=ext_nr, extf = 1.1)

# column names for both background_points
colnames(bck) = c('lon', 'lat') # name columns for background points
colnames(nr_bck) = c('lon', 'lat')

# create sp layer with coordinate reference system
bck <- SpatialPoints(bck, crs(slp))
nr_bck <- SpatialPoints(nr_bck, crs(slp))

plot(bck, pch = 20) # check background points
plot(nr_bck, pch = 20, col = 'red', add = TRUE)
# -------------- Partition Data --------------
# group background points/pseduo-absences - for global range
group = kfold(bck, k=5)
# group pseudo-absences for training model performances
bck_test = bck[group== 1,]
bck_train = bck[group!= 1,]

# group background points/pseduo-absences - for native range
nr_group = kfold(nr_bck, k=5)
# group pseudo-absences for training model performances
nr_bck_test = nr_bck[nr_group== 1,]
nr_bck_train = nr_bck[nr_group!= 1,]

# ----------------------- MAXENT MODEL --------------------------------
#------------------------- FIT MAXENT -----------------------#
# fit maxent - using stacked raster variables (one global and one cropped to native range (max model nr))
maxent_model_wrld <- maxent(predictor, atulip_train_wrld) # global model with global distribution
maxent_model <- maxent(predictor, atulip_train) # global predictor and native range presence
max_model_nr <- maxent(nr_predictor, atulip_train) # native range predictor and presence

# train maxent model, use absence points, and biome as factor and one restricted to NR
prob_model_wrld <- maxent(predictor, atulip_train_wrld, bck_train, factor = "biome.tif") # global model with global pres/abs
prob_model <- maxent(predictor, atulip_train, nr_bck_train, factor = "biome.tif") # global pred with nr pres/abs
prob_model_nr <- maxent(nr_predictor, atulip_train, nr_bck_train, factor = "biome.tif") # nr pred, pres/abs

# plot each variable
par(mfrow = c(2,3))
plot(maxent_model_wrld, main = "Maxent Global Model - Presence Only")
plot(maxent_model, main = "Maxent Model NR Pres - Presence Only")
plot(max_model_nr, main = "Maxent Model NR - Presence Only")

plot(prob_model_wrld, main = "Maxent Global Model - Pres/Abs")
plot(prob_model, main = "Maxent Model Global Pred - Pres/Abs")
plot(prob_model_nr, main = "Maxent Model NR - Pres/Abs")

title <- names(predictor)
ls <- 1:9

par(mfrow = c(2,2), mar = rep(2, 4))

# run response curve - response curves are so different
for (i in ls){
  x <- title[i]
  response(maxent_model, var = x, main = x)
}

for (i in ls){
  x <- title[i]
  response(max_model_nr, var = x, main = x)
}

for (i in ls){
  x <- title[i]
  response(prob_model, var = x, main = x)
}

for (i in ls){
  x <- title[i]
  response(prob_model_nr, var = x, main = x)
}

par(mfrow = c(1,1))
# ---------- MODEL EVALUATION -------------#
# PRESENCE ONLY TRAINED MODEL
# run evaluation - using predictor variables 
# global prediction (on presence only)
e_wrld = evaluate(atulip_test_wrld, bck_test, maxent_model_wrld, predictor)
e_wrld

e = evaluate(atulip_test, nr_bck_test, maxent_model, predictor)
e

e_suit_nr = evaluate(atulip_test, nr_bck_test, max_model_nr, nr_predictor)
e_suit_nr

# PRES/ABS TRAINED MODEL (nr pres/abs only, all NR data)
# evaluate probability model based on absence data
e_prob_wrld = evaluate(atulip_test_wrld, bck_test, prob_model_wrld, predictor)
e_prob_wrld

#evaluate probability model pres/abs within native range and global predictor
e_prob = evaluate(atulip_test, nr_bck_test, prob_model, predictor)
e_prob

# evaluate probability model all predictors, pres/abs within native range
e_prob_nr = evaluate(atulip_test, nr_bck_test, prob_model_nr, nr_predictor)
e_prob_nr

par(mfrow = c(1,2))
# plot ROC - and AUC
plot(e, 'ROC', main='ROC Curve - Maxent Pres Only')
plot(e, 'TPR', main='TPR Curve - Maxent Pres Only')

plot(e_prob_nr, 'ROC', main='ROC Curve - Maxent Pres Only NRP')
plot(e_prob_nr, 'TPR', main='TPR Curve - Maxent Pres Only NRP')


plot(e_prob, 'ROC', main='ROC Curve - Maxent Pres/Abs GP')
plot(e_prob, 'TPR', main='TPR Curve - Maxent Pres/Abs GP')

plot(e_nr, 'ROC', main='ROC Curve - Maxent Pres/Abs NRP')
plot(e_nr, 'TPR', main='TPR Curve - Maxent Pres/Abs NRP')

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
pic_predictor <- projectRaster(predictor, crs = pdc) # fairly long processing (1 day)

# crop reprojected rasters based on extent of pacific region shp
pic_predictor <- crop(pic_predictor, pic_ext)

# convert global atulip test/train data for evaluation
atulip_pic_train <- spTransform(atulip_train_wrld, crs(pdc))
atulip_pic_test <- spTransform(atulip_test_wrld, crs(pdc))

# crop reprojected data to pacific extent to run evaluation
atulip_pic_train <- crop(atulip_pic_train, pic_ext)
atulip_pic_test <- crop(atulip_pic_test, pic_ext)


# ---------------------------- RUN PREDICTION ON PACIFIC REGION -----------------------------
# Presence Only Models
#predict the entire training dataset
atulip_pred_wrld <- predict(maxent_model_wrld, pic_predictor, ext = pic_ext)

atulip_pred <- predict(maxent_model, pic_predictor, ext=pic_ext)

atulip_pred_nr <- predict(max_model_nr, pic_predictor, ext=pic_ext)

# Presence and Absence Fitted Models
# predict the entire dataset to Pacific region using nr data 
atulip_prob_wrld <- predict(prob_model_wrld, pic_predictor, ext=pic_ext)

atulip_prob <- predict(prob_model, pic_predictor, ext=pic_ext)

# run nr model (native range occurrence fit to native range and predicted to pacific extent)
atulip_nr_prob <- predict(prob_model_nr, pic_predictor, ext = pic_ext)

# setting color ramp
clr <- RColorBrewer::brewer.pal(9, "YlOrRd")
rdcl <- RColorBrewer::brewer.pal(9, 'Reds')
crmp <- colorRampPalette(clr)
rdcrmp <- colorRampPalette(rdcl)
par(mfrow=c(1,1))

# plot maps

par(mfrow=c(2,3))
plot(atulip_pred_wrld, col=rdcrmp(10), main = "Pres Only - 30s African Tulip Global Suitability for Pacific Region")
plot(atulip_pred, col=rdcrmp(10), main = "Pres Only GPred - 30s African Tulip Suitability for Pacific Region")
plot(atulip_pred_nr, col=rdcrmp(10), main = "Pres Only NRP - 30s African Tulip Suitability for Pacific Region")

plot(atulip_prob_wrld, col=rdcrmp(10), main = "Pres/Abs - 30s African Tulip Global Probabilty for Pacific Region")
plot(atulip_prob, col=rdcrmp(10), main = "Pres/Abs GPred - 30s Africant Tulip Probability for Pacific Region")
plot(atulip_nr_prob, col=rdcrmp(10), main = "Pres/Abs NRP - 30s African Tulip Probability for Pacific Region")


#### --------------------- BOYCE INDEX -------------###
# run boyce index for model evaluation
# remotes::install_github('adamlilith/enmSdm', dependencies=TRUE)
library(enmSdm)

# PRESENCE ABSENCE MODELS
# Pres/Abs NRP - African Tulip Probability for Pacific Region
# (native range occurrence fit to native range and predicted to pacific extent)
atulip_nr_prob_1 <- extract(atulip_nr_prob, atulip_pic_test)

# crop testing background absence points
atulip_bck_test <- spTransform(bck_test, crs(pdc))
atulip_bck_test <- crop(atulip_bck_test, pic_ext)
atulip_nr_prob_abs <- extract(atulip_nr_prob, atulip_bck_test)

contBoyce(pres = atulip_nr_prob_1, 
          contrast = atulip_nr_prob_abs,
          na.rm = T)

# Pres/Abs NRP - African Tulip Probability for Pacific Region
atulip_prob_1 <- extract(atulip_prob, atulip_pic_test)
atulip_prob_abs <- extract(atulip_prob, atulip_bck_test)

contBoyce(pres = atulip_prob_1, 
          contrast = atulip_prob_abs,
          na.rm = T)

# Global Prob Model 
atulip_prob_wrld_1 <- extract(atulip_prob_wrld, atulip_pic_test)
atulip_prob_wrld_abs <- extract(atulip_prob_wrld, atulip_bck_test)

contBoyce(pres = atulip_prob_wrld_1, 
          contrast = atulip_prob_wrld_abs,
          na.rm = T)

# --------------PRESENCE ONLY MODELS
# global predictors and presence points
atulip_pred_wrld_1 <- extract(atulip_pred_wrld, atulip_pic_test) # predicted presence values
atulip_pred_wrld_abs <- extract(atulip_pred_wrld, atulip_bck_test) # predicted absence values

contBoyce(pres = atulip_pred_wrld_1,
          contrast = atulip_pred_wrld_abs,
          na.rm = T)

# global predictors and native range presence points
atulip_pred_1 <- extract(atulip_pred, atulip_pic_test)
atulip_pred_abs <- extract(atulip_pred, atulip_bck_test)

contBoyce(pres = atulip_pred_1,
          contrast = atulip_pred_abs,
          na.rm = T)


# native predictors and native range presence points
atulip_pred_nr_1 <- extract(atulip_pred_nr, atulip_pic_test)
atulip_pred_nr_abs <- extract(atulip_pred_nr, atulip_bck_test)

contBoyce(pres = atulip_pred_nr_1,
          contrast = atulip_pred_nr_abs,
          na.rm = T)

# ------------------- Plot higher performing models to country -----------------
# country level - Fiji
install.packages('rnaturalearthhires', repos = "http://packages.ropensci.org", type = "source")
require(rnaturalearth)
library(rnaturalearthhires)
fj.p <- ne_countries(country = 'fiji', scale = 'large')
fj.p <- spTransform(fj.p, crs(pdc))
plot(fj.p)

fj <- extent(fj.p)
plot(fj)

par(mfrow=c(2,2))
# pres only global model
plot(atulip_pred_wrld, col=rdcrmp(10), ext = fj, main = "Pres Only - African Tulip Global Suitability for Fiji")
plot(fj.p, add = TRUE)
# native range presence only (native range predictors)
plot(atulip_pred_nr, col=rdcrmp(10), ext = fj, main = "Pres Only NRP Probabilty for Fiji")
plot(fj.p, add = TRUE)
# pres/abs global model
plot(atulip_prob_wrld, col=rdcrmp(10), ext = fj, main = "Pres/Abs - African Tulip Global Probabilty for Fiji")
plot(fj.p, add = TRUE)
# pres/abs native range model 
plot(atulip_nr_prob, col=rdcrmp(10), ext = fj, main = "Pres/Abs NRP - African Tulip Probability for Fiji")
plot(fj.p, add = TRUE)
# -------- SAMOA -------------
# country level - Samoa
ws.p <- ne_countries(country = 'samoa', scale = 'large')
ws.p <- spTransform(ws.p, crs(pdc))
plot(ws.p)
# get extent of Samoa
ws <- extent(ws.p)
plot(ws)

par(mfrow=c(2,2))
# pres only global model
plot(atulip_pred_wrld, col=rdcrmp(10), ext = ws, main = "Pres Only - African Tulip Global Suitability for Samoa")
plot(ws.p, add = TRUE)
# native range presence only (native range predictors)
plot(atulip_pred_nr, col=rdcrmp(10), ext = ws, main = "Pres Only NRP Probabilty for Samoa")
plot(ws.p, add = TRUE)
# pres/abs global model
plot(atulip_prob_wrld, col=rdcrmp(10), ext = ws, main = "Pres/Abs - African Tulip Global Probabilty for Samoa")
plot(ws.p, add = TRUE)
# pres/abs native range model 
plot(atulip_nr_prob, col=rdcrmp(10), ext = ws, main = "Pres/Abs NRP - African Tulip Probability for Samoa")
plot(ws.p, add = TRUE)


# -------- TONGA -------------
# country level - Tonga
to.p <- ne_countries(country = 'tonga', scale = 'large')
to.p <- spTransform(to.p, crs(pdc))
plot(to.p)
# plot extent of Tonga
to <- extent(to.p)
plot(to)

par(mfrow=c(2,2))
# pres only global model
plot(atulip_pred_wrld, col=rdcrmp(10), ext = to, main = "Pres Only - African Tulip Global Suitability for Tonga")
plot(to.p, add = TRUE)
# native range presence only (native range predictors)
plot(atulip_pred_nr, col=rdcrmp(10), ext = to, main = "Pres Only NRP Probabilty for Tonga")
plot(to.p, add = TRUE)
# pres/abs global model
plot(atulip_prob_wrld, col=rdcrmp(10), ext = to, main = "Pres/Abs - African Tulip Global Probabilty for Tonga")
plot(to.p, add = TRUE)
# pres/abs native range model 
plot(atulip_nr_prob, col=rdcrmp(10), ext = to, main = "Pres/Abs NRP - African Tulip Probability for Tonga")
plot(to.p, add = TRUE)


