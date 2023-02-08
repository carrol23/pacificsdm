# ------ Set Working Directory --------
setwd('C:/Data/niche-modelling/SDM_22/merremia/30sec/data')

rm(list=ls()) # clean workspace

# install and load packages
install.packages(c("dismo","maptools","maps", "rJava", "tidyverse", "ggplot2", "raster","sp"))
library(dismo)
library(maptools) # retires in 2023
library(maps)
library(ggplot2)

library(sp)
library(raster)
library(rJava)

library(rgbif)
library(terra)
library(tidyverse)
library(CoordinateCleaner)

# ---------- DOWNLOAD SPECIES DATA -------
# run occurence search using altnerative name (combined data under Merremia peltata scientific name)
merr_alt = occ_search(scientificName = "Decalobanthus peltatus", limit = 10000, hasCoordinate = TRUE) # 691 observations

# select data
merr_alt <- merr_alt$data
# add in ALA dataset - australian natural history museum
merr_csv <- read.csv('./data/ala/records-2022-11-20.csv')

# view data for both imported
View(merr_alt)
View(merr_csv)

# scale data down - remove columns not useful (only select the following)
merr_alt <- merr_alt %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)

# remove records without coordinates - records should be the same
merr_alt <- merr_alt %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# identify duplicates based on lat and long coordinates
merr_dup = duplicated(merr_alt[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 256 duplicates (out of 691 observations)
length(merr_dup[merr_dup == TRUE])

# remove duplicates - returns 435 observations
merr_data <- merr_alt[!merr_dup,]

# remove records without coordinates and duplicates for csv loaded
merr_csv <- merr_csv %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

merr_csv_dup <- duplicated(merr_csv[,c("decimalLongitude", "decimalLatitude")])
length(merr_csv_dup[merr_csv_dup == TRUE]) # get number of duplicated = 98
merr_csv <- merr_csv[!merr_csv_dup,]

# scale data down - remove columns not useful (only select the following)
merr_csv <- merr_csv %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         dataResourceUid, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)


# plot remaining occurrences
wrld = borders("world", colour = "gray50", fill = "gray80")

mplot <- ggplot() +  wrld +
  geom_point(data = merr_data, aes(x=decimalLongitude, y=decimalLatitude),
             colour = "yellow", size = 3, alpha = 0.1) +
  geom_point(data = merr_csv, aes(x = decimalLongitude, y = decimalLatitude), colour = "blue", size = 1.5) + 
  coord_sf()
  

mplot

# ---------------------- COORDINATE CLEANER TESTS ---------------------
# use coordinate cleaner to flag problem records
library(countrycode)

# convert country code ISO2 to ISO3
merr_data$countryCode <- countrycode(merr_data$countryCode, 
                                       origin = 'iso2c',
                                       destination = 'iso3c')
merr_csv$countryCode <- countrycode(merr_data$countryCode, 
                                       origin = 'iso2c',
                                       destination = 'iso3c')

# flag coordinates not 'clean' based on tests
flags <- clean_coordinates(x = merr_data,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'outliers', 'seas', 'countries'))

# flag coordinates not 'clean' based on tests
flags_csv <- clean_coordinates(x = merr_csv,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'seas')) #countries


summary(flags) # 180 flagged records
summary(flags_csv) # 73 flagged records
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
plot(flags_csv, lon = "decimalLongitude", lat = "decimalLatitude")

# exclude flagged records
merr_occ <- merr_data[flags$.summary,] # 262
merr_occ_csv <- merr_csv[flags_csv$.summary,] #565

# merge datasets
merr_occ = subset(merr_occ, select=-c(gbifID))
merr_occ_csv = subset(merr_occ_csv, select = -c(dataResourceUid))

merr_occ <- rbind(merr_occ, merr_occ_csv)

# identify duplicates based on lat and long coordinates
merr_occ_dup = duplicated(merr_occ[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 3 duplicates (out of 827 observations)
length(merr_occ_dup[merr_occ_dup == TRUE])

# remove duplicates - returns 435 observations
merr_occ <- merr_occ[!merr_occ_dup,]

# ----- METADATA TEST -------
# remove records with low precision
hist(merr_occ$coordinateUncertaintyInMeters / 1000, breaks = 20)

# removed coordinates with uncertainty up to 100km
merr_occ <- merr_occ %>%
  filter(coordinateUncertaintyInMeters/1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# identify basis of records - remove fossil, living & preserved specimen 
table(merr_occ$basisOfRecord)

merr_data_filt <- filter(merr_occ, basisOfRecord == "HUMAN_OBSERVATION" |
                        basisOfRecord == "OBSERVATION" |
                        basisOfRecord == "PRESERVED_SPECIMEN" |
                        basisOfRecord == "OCCURENCE")

table(merr_data_filt$individualCount)

# remove individual counts that equate to 0 and higher than 100
merr_data_filt <-merr_data_filt %>%
  filter(individualCount > 0 | is.na(individualCount)) %>%
  filter(individualCount < 99 | is.na(individualCount))

# exclude old records- 1925 remove
table(merr_data_filt$year)

merr_occ <- merr_data_filt %>%
  filter(year > 1945)

# check family name
table(merr_occ$family)

# check taxon rank
table(merr_occ$taxonRank)

# Split occurence data within defined Native Range
library(rgdal)
nr <- readOGR("C:/Data/niche-modelling/SDM_22/merremia/30sec/nr_merr.shp")

# convert to df and add in columns (id and species)
nr_df <- as.data.frame(nr, row.names = T)
nr_df$id <- c("A")
nr_df$species <- c("Decalobanthus peltatus")

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
  x = merr_occ,
  range = nr_shp,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged"
)

# devil data within the range
data_fin <- merr_occ[range_flags,]

range +
  geom_point(data = merr_occ, aes(x=decimalLongitude, y = decimalLatitude), color = 'blue') +
  geom_point(data = data_fin, aes(x=decimalLongitude, y=decimalLatitude), color = 'red')

# ------------- Problematic datasets -------------
# ddmm to dd.dd conversion error - 0 records are flagged
out.ddmm <- cd_ddmm(merr_occ, lon = "decimalLongitude",
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

# ----- DROP BIOCLIMATIC VARIABLES -----------------
# Limit set of bioclimatic predictors
# available list of bioclimatic indicators - https://www.worldclim.org/data/bioclim.html 
currentEnv = dropLayer(currentEnv.ras, c("bio2", "bio3", "bio4", "bio10","bio11", 
                                     "bio13", "bio14", "bio15", "bio18", "bio19"))

# write each raster to file
writeRaster(stack(currentEnv), names(currentEnv), bylayer = TRUE,
            format = 'GTiff', overwrite=TRUE)

# ------- ELEVATION VARIABLES ---------- #
# obtain other environmental variables - elevation, slope, hillshade
# download elevation variables
elevation_wrld <- worldclim_global(var="elev", res=0.5, path = ".") # alt from srtm
#convert to raster
elevation_wrld <- as(elevation_wrld, "Raster")
raster::writeRaster(elevation_wrld, filename='world_elevation.tif', format='GTiff', overwrite=TRUE)

# run topographic analysis
# calculate slope
slp = terra::terrain(elevation_wrld, v = 'slope', unit = 'radians', neighbors = 8, 
                     filename = 'slope_wrld.tif', overwrite = TRUE)

# convert spat raster to raster layer
slp <- as(slp, "Raster")
plot(slp)

# calculate aspect
aspect = terrain(elevation_wrld, opt='aspect', filename='world_aspect.tif', overwrite=TRUE)

# -------- Categorical Variable - Terrestrial Ecoregions --------
# read in raster biome
biome <- raster("C:/Data/niche-modelling/SDM_22/merremia/data/biome/wrld_biome.tif")
# clip biome to extent - remove polar regions 

biome <- crop(stck_biome, ext, mask = T)
# resample biome layer to match extent, crs and resolution of elevation_wrld
stck_biome <- resample(biome, elevation_wrld, method = "ngb")

# assign projection from slp layer to biome dataset
writeRaster(biome, filename="biome.tif", format='GTiff', overwrite=TRUE)

# --------MERGE DOWNLOADED VARIABLES----------
# stack raster layers
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

# --------------------Check CORRELATION ------------------
library(sdmpredictors)

# run correlation co-efficient
p_mat <- pearson_correlation_matrix(stck, same_mask = FALSE)
View(p_mat) # remove layers with high correlation > 0.8 (bio 12 and bio 16)

# drop highly coorelated predicts
stck = dropLayer(stck, c("bio1.tif", "bio6.tif", "bio17.tif", "world_aspect.tif", 
                         "preciptn.tif", "tempavg.tif", "tempmin.tif", "bio13.tif"))


# ------- GENERATE NATIVE RANGE PREDICTOR --------------------
# generate extent from native range
ext_nr <- raster::extent(nr_shp)
nr_predictor <- crop(stck, ext_nr, snap = 'near', datatype = NULL) # outputs raster brick
names(nr_predictor)

plot(nr_predictor)

# ----- PARTITION DATA ---------
# SPLIT PRESENCE POINTS
# 20% sample for testing
merr_df = cbind.data.frame(merr_occ$decimalLongitude, merr_occ$decimalLatitude)
merr_df <- SpatialPoints(merr_df, crs(slp))
fold <- kfold(merr_df, k=5)
merr_test <- merr_df[fold == 1,]
merr_train <- merr_df[fold!=1,]

# native range data 
nr_data <- data_fin[,2:3]
nr_data <- SpatialPoints(nr_data, crs(slp))
fold <- kfold(nr_data, k = 5)
nr_test <- nr_data[fold == 1,] # 20% of total occ data from native range
nr_train <- nr_data[fold != 1,] # 80 of total occ data from native range

# ------------------ABSENCE POINTS------------------- 
# generate background points/pseudo-abs
merr_bck = randomPoints(stck, n = 300,p = merr_df, ext = ext, extf = 1.1)
nr_bck = randomPoints(nr_predictor, n = 200, p = nr_data, ext = nr, extf = 1.1)

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

#---------- FIT MAXENT ------------#
# fit maxent - using stacked raster predictor variables
# PRES ONLY MODELS
model <- maxent(stck, merr_train)
model_nrp <- maxent(stck, nr_train) 
nr_model <- maxent(nr_predictor, nr_train)

# PRES/ABS MODELS
prob_model <- maxent(stck, merr_train, merr_bck_train, factor = "biome.tif")
prob_model_nrp <- maxent(stck, nr_train, merr_nr_bck_train, factor = "biome.tif")
nr_prob_model <- maxent(nr_predictor, nr_train, merr_nr_bck_train, factor = "biome.tif")

par(mfrow = c(1,3))
plot(model, main = "Presence Only Global Merremia Maxent Model")
plot(model_nrp, main = "Pres Only Native Range Occurence Model")
plot(nr_model, main = "Pres Only Native Range Restricted Model")

par(mfrow = c(1,3))
plot(prob_model, main = "Pres-Abs Global Merremia Maxent Model")
plot(prob_model_nrp, main = "Pres-Abs Native Range Occurence Model")
plot(nr_prob_model, main = "Pres-Abs Native Range Restricted Model")

# run response curve
response(model)
response(model_nrp)
response(nr_model)

response(prob_model)
response(prob_model_nrp)
response(nr_prob_model)

# remove low perforing variables - re-run model model
# drop highly coorelated predicts
stck = dropLayer(stck, c("bio16.tif", "bio9.tif", "bio8.tif"))
nr_predictor = dropLayer(nr_predictor, c("bio16.tif", "bio9.tif", "bio8.tif"))

# ---------- EVALUATE TRAINED MODEL -------------#
# run evaluation - using predictor variables
e1 = evaluate(merr_test, merr_bck_test, model, stck) # low number of test points
e2 = evaluate(nr_test, merr_nr_bck_test, model_nrp, stck) # 
e3 = evaluate(nr_test, merr_nr_bck_test, nr_model, nr_predictor) # 

e1
e2
e3

prob_e1 = evaluate(merr_test, merr_bck_test, prob_model, stck)
prob_e2 = evaluate(nr_test, merr_nr_bck_test, prob_model_nrp, stck)
prob_e3 = evaluate(nr_test, merr_nr_bck_test, nr_prob_model, nr_predictor)

prob_e1
prob_e2
prob_e3

# plot ROC - and AUC
par(mfrow = c(1,2))
plot(e1, 'ROC', main = 'ROC Curve - Maxent Current Predictor Variables')
plot(e1, 'TPR', main = 'TPR Curve - Maxent Current Preditor Variables')

plot(prob_e1, 'ROC', main = 'ROC Curve - Maxent Current Predictor Variables')
plot(prob_e1, 'TPR', main = 'TPR Curve - Maxent Current Preditor Variables')

# ------------------ RUN PREDICTION ON LARGER AREA ---------------- 
# --------------------- CROP DATA TO PACIFIC REGION ------------
# restrict data to pacific/oceania region (def as region 5)
data(wrld_simpl)
pc <- subset(wrld_simpl, wrld_simpl$REGION == 9)
# reproject pc first
pdc <- crs.latlong <- crs("+init=epsg:3832")
pc_rpj <- spTransform(pc, crs(pdc))
plot(pc_rpj)

# get extent
pic_ext <- extent(pc_rpj)

# ------------------- CROP PREDICTOR VARIABLES ------------------
# get predictor variables for PIC region
#pic_predictor <- projectRaster(stck, crs = pdc) # long processing 3 days (undertaken in QGIS as quicker alternative)

# get stck reprojected in qgis
setwd('C:/Data/niche-modelling/SDM_22/merremia/30sec/pic_crop')
pic_datafiles = Sys.glob("*.tif")
pic_predictor = stack()
for (i in 1:NROW(pic_datafiles)){
  tempraster = raster(pic_datafiles[i]) #iterate each file and store as temporary raster
  pic_predictor = stack(pic_predictor, tempraster)
}

# get list of names from file
r.ls <- (list.files(pattern = "tif$"))

# rename stack with layer names
names(pic_predictor) <- r.ls
# check names
names(pic_predictor)

# crop reprojected rasters based on extent of pacific region shp
pic_predictor <- crop(pic_predictor, pic_ext)

# --------------------- RUN PREDICTION -------------------------------
# predict  trained models - pres only
merr_predict <- predict(model, pic_predictor, ext = pic_ext)
merr_predict_nrp <- predict(model_nrp, pic_predictor, ext = pic_ext)
merr_predict_nr <- predict(nr_model, pic_predictor, ext = pic_ext)

# predict trained model - pres/abs 
merr_prob_predict <- predict(prob_model, pic_predictor, ext = pic_ext)
merr_prob_predict_nrp <- predict(prob_model_nrp, pic_predictor, ext = pic_ext)
merr_prob_predict_nr <- predict(nr_prob_model, pic_predictor, ext = pic_ext)

#--------------------- COLOR RAMP ------------------
# setting color ramp
clr <- RColorBrewer::brewer.pal(9, "YlOrRd")
rdcl <- RColorBrewer::brewer.pal(9, 'Reds')
crmp <- colorRampPalette(clr)
rdcrmp <- colorRampPalette(rdcl)

par(mfrow=c(2,3))
# plot predicted model - presence only models
plot(merr_predict, col = rdcrmp(10), main = 'Global Predicted Suitability for Merremia - Pres Only')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_predict, 'merremia_predicted_global.tif')

plot(merr_predict_nrp, col = rdcrmp(10), main = 'NRP Predicted Suitability for Merremia - Pres Only')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_predict_nrp, 'merremia_predicted_nrp.tif')

plot(merr_predict_nr, col = rdcrmp(10), main = 'Native Range Predicted Suitability for Merremia - Pres Only')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_predict_nrp, 'merremia_predicted_nr.tif')

# plot predicted model - pres/abs models
plot(merr_prob_predict, col = rdcrmp(10), main = 'Predicted Probability for Merremia - Pres/Abs')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_prob_predict, 'merremia_prob_global.tif')

plot(merr_prob_predict_nrp, col = rdcrmp(10), main = 'NRP Predicted Probability for Merremia - Pres/Abs')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_prob_predict_nrp, 'merremia_prob_nrp.tif')

plot(merr_prob_predict_nr, col = rdcrmp(10), main = 'Native Range Predicted Probability for Merremia - Pres/Abs')
maps::map(database = 'world', fill = FALSE, add = TRUE)
writeRaster(merr_prob_predict_nr, 'merremia_prob_global_nr.tif')

# ----------------------------- BOYCE INDEX --------------------------#
# run boyce index for model evaluation
# remotes::install_github('adamlilith/enmSdm', dependencies=TRUE)
library(enmSdm)

# CREATE PRESENCE POINT DATA
merr_df_pic <- spTransform(merr_df, crs(pdc))
merr_df_pic <- crop(merr_df_pic, pic_ext) # results in 133 points within the pacific region

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
contBoyce(pres = merr_prob_pres, 
          contrast = merr_prob_abs,
          na.rm = T)

# Native range presence and global predictors
merr_prob_pres <- raster::extract(merr_predict_nrp, merr_df_pic)
merr_prob_abs <- raster::extract(merr_predict_nrp, merr_bck_pic)
contBoyce(pres = merr_prob_pres, 
          contrast = merr_prob_abs,
          na.rm = T)

# native range presence and predictor
merr_prob_pres <- raster::extract(merr_predict_nr, merr_df_pic)
merr_prob_abs <- raster::extract(merr_predict_nr, merr_bck_pic)
contBoyce(pres = merr_prob_pres, 
          contrast = merr_prob_abs,
          na.rm = T)

# PRES ABSENCE MODEL PREDICTION
# Global
merr_prob_1 <- raster::extract(merr_prob_predict, merr_df_pic)
merr_prob_abs1 <- raster::extract(merr_prob_predict, merr_bck_pic)
contBoyce(pres = merr_prob_1, 
          contrast = merr_prob_abs1,
          na.rm = T)
# Native range presence only
merr_prob_1 <- raster::extract(merr_prob_predict_nrp, merr_df_pic)
merr_prob_abs1 <- raster::extract(merr_prob_predict_nrp, merr_bck_pic)
contBoyce(pres = merr_prob_1, 
          contrast = merr_prob_abs1,
          na.rm = T)
# Native range presence and predictors
merr_prob_1 <- raster::extract(merr_prob_predict_nr, merr_df_pic)
merr_prob_abs1 <- raster::extract(merr_prob_predict_nr, merr_bck_pic)
contBoyce(pres = merr_prob_1, 
          contrast = merr_prob_abs1,
          na.rm = T)

# ------------------- Plot higher performing models to country -----------------
# country level - Fiji
fj.p <- subset(pc_rpj, pc_rpj$ISO2 == 'FJ')
fj <- extent(fj.p)

par(mfrow=c(2,3))
# pres only global model
plot(merr_predict, col=rdcrmp(10), ext = fj, main = "Pres Only - Merremia Suitability for Fiji")
plot(merr_predict_nrp, col=rdcrmp(10), ext = fj, main = "Pres Only - NR Presence - Merremia Suitability for Fiji")
plot(merr_predict_nr, col=rdcrmp(10), ext = fj, main = "Pres Only - NR - Merremia Suitability for Fiji")

# pres/abs only model
plot(merr_prob_predict, col=rdcrmp(10), ext = fj, main = "Pres/Abs - Merremia Probability for Fiji")
plot(merr_prob_predict_nrp, col=rdcrmp(10), ext = fj, main = "Pres/Abs - NR Presence - Merremia Probability for Fiji")
plot(merr_prob_predict_nr, col=rdcrmp(10), ext = fj, main = "Pres/Abs - NR - Merremia Probability for Fiji")


# -------- SAMOA -------------
install.packages('rnaturalearthhires', repos = "http://packages.ropensci.org", type = "source")
require(rnaturalearth)
library(rnaturalearthhires)
# country level - Samoa
ws.p <- ne_countries(country = 'samoa', scale = 'large')
ws.p <- spTransform(ws.p, crs(pdc))
plot(ws.p)
# get extent of Samoa
ws <- extent(ws.p)
plot(ws)

par(mfrow=c(2,3))
# pres only global model
plot(merr_predict, col=rdcrmp(10), ext = ws, main = "Pres Only - Merremia Suitability for Samoa")
plot(merr_predict_nrp, col=rdcrmp(10), ext = ws, main = "Pres Only - NR Presence - Merremia Suitability for Samoa")
plot(merr_predict_nr, col=rdcrmp(10), ext = ws, main = "Pres Only - NR - Merremia Suitability for Samoa")

# pres/abs only model
plot(merr_prob_predict, col=rdcrmp(10), ext = ws, main = "Pres/Abs - Merremia Probability for Samoa")
plot(merr_prob_predict_nrp, col=rdcrmp(10), ext = ws, main = "Pres/Abs - NR Presence - Merremia Probability for Samoa")
plot(merr_prob_predict_nr, col=rdcrmp(10), ext = ws, main = "Pres/Abs - NR - Merremia Probability for Samoa")


# -------- TONGA -------------
# country level - Tonga
to.p <- ne_countries(country = 'tonga', scale = 'large')
to.p <- spTransform(to.p, crs(pdc))
plot(to.p)
# plot extent of Tonga
to <- extent(to.p)
plot(to)

par(mfrow=c(2,3))
# pres only global model
plot(merr_predict, col=rdcrmp(10), ext = to, main = "Pres Only - Merremia Suitability for Tonga")
plot(merr_predict_nrp, col=rdcrmp(10), ext = to, main = "Pres Only - NR Presence - Merremia Suitability for Tonga")
plot(merr_predict_nr, col=rdcrmp(10), ext = to, main = "Pres Only - NR - Merremia Suitability for Tonga")

# pres/abs only model
plot(merr_prob_predict, col=rdcrmp(10), ext = to, main = "Pres/Abs - Merremia Probability for Tonga")
plot(merr_prob_predict_nrp, col=rdcrmp(10), ext = to, main = "Pres/Abs - NR Presence - Merremia Probability for Tonga")
plot(merr_prob_predict_nr, col=rdcrmp(10), ext = to, main = "Pres/Abs - NR - Merremia Probability for Tonga")

