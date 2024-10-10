# PRE-PROCESS DATA FOR MODEL
# Load in occurence data and perform data thinning
# ------ Set Working Directory --------
setwd('C:/Data/sdm/C2_SDM/')
rm(list = ls()) # clean workspace

library(sf)
library(raster)
library(terra)
library(fuzzySim)
library(devtools)
# devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")
library(mecofun)

# ----------- LOAD OCCURRENCE DATA -----------
# read in occurrence points that have been cleaned and exported to local folder
atulip_occ = read.csv("C:/Data/sdm/C2_SDM/occurence/at/comb_at.csv", sep = ',')

# convert presence points to spatial points
occ_df = cbind.data.frame(atulip_occ$decimalLongitude, atulip_occ$decimalLatitude)
names(occ_df)[names(occ_df) == 'occ$decimalLongitude'] <- 'lon'
names(occ_df)[names(occ_df) == 'occ$decimalLatitude'] <- 'lat'
occ_sp = SpatialPoints(occ_df)
crs(occ_sp) <- crs("+init=epsg:4326")

# ------------ GET PSEUDO-ABSENCE DATA -----------
# generate background points as pseudo-absence points
# load in currentEnv to get extent
currentEnv <- stack('./clim/cor_stack.tif')
names(currentEnv) <- c("bio1", "bio2", "bio3", "bio7",
                           "bio12", "bio14", "bio18", "elevation")

ext = extent(currentEnv) #get extent
abs = randomPoints(currentEnv, n = 3000,p = occ_sp, ext = ext, extf = 1.0)

# column names for both background_points
colnames(abs) = c('lon', 'lat')
# create sp layer with coordinate reference system for global generated bckgrd points
abs_sp <- SpatialPoints(abs)
crs(abs_sp) <- crs("+init=epsg:4326")

# PLOT OCCURRENCE AND ABSENCE
plot(occ_df, pch = 20, col = 'blue')
plot(abs, pch = 20, col = 'red', add = TRUE)

# THIN ALL DATA
# thin occurrence data based on the pixel resolution of raster layer (stck)
rec_df <- gridRecords(currentEnv, occ_sp, abs_sp, na.rm = TRUE)
rec_sp_xy <- SpatialPointsDataFrame(coords = rec_df[,c(2,3)], data = rec_df, proj4string = crs(currentEnv))
rec_sp_xy <- rec_sp_xy[,c(1,5:12)] # remove cells column (not needed)

write.csv(cbind(coordinates(rec_sp_xy), rec_sp_xy@data),
          file = "C:/Data/sdm/C2_SDM/occurence/at/at_thinned_comb_all.csv",
          row.names = FALSE )

# ------- CROP OCCURENCE DATA TO PACIFIC REGION ------
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

# get extent
pic_ext <- extent(pc_rpj)

# GET PACIFIC OCC AND ABS DATA
# reproject the occurrence/abs data and crop to Pacific
pic_sp <- spTransform(rec_sp_xy, pdc)
pic_sp <- crop(pic_sp, pic_ext)
# returns 646 points (pres/abs)

write.csv(cbind(coordinates(pic_sp), pic_sp@data),
          file = "C:/Data/sdm/C2_SDM/occurence/at/at_thinned_comb_pic.csv",
          row.names = FALSE )