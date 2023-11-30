#  ---------------------SDM ENVIRONMETNAL PREDICTORS DATA  -------------------------
# 1 download data from worldclim, biome and soil data
# 2 clean and merge datasets
# 3 run pearsons correlation to drop variables
# 4 write predictors to folder

# ------ Set Working Directory --------
setwd('C:/Data/niche-modelling/SDM_22/merremia/30sec/data')

rm(list=ls()) # clean workspace

# ---------------- INSTALL PACKAGES -------------
install.packages(c('dismo', 'sp', 'maps', 'tidyverse', 'ggplot2', 'geodata'))
library(dismo)
library(maps)
library(ggplot2)
library(sp)
library(tidyverse)

# ------- Access and Download Environmental Data ----------
# -------------- WORLDCLIM DATA -------------
# get WORLDCLIM data
library(geodata)
currentEnv = worldclim_global(var = "bio", res = 0.5, path = ".") # download 0.5 mins (30sec data0)

# convert spatraster to raster
currentEnv.ras <- as(currentEnv, "Raster")

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
            format = 'GTiff', overwrite = TRUE)

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

# ------------SOIL DATA ----------------
# read in soil data - GSOC soil
gsoc <- raster("C:/Data/niche-modelling/global/GSOCmap/GSOC.tif")
gsoc <- crop(gsoc, extent(slp), mask  = T)
soil <- raster::projectRaster(gsoc, slp, method = "bilinear")
writeRaster(soil, filename ="gsoc.tif", format ='GTiff', overwrite = TRUE)

# -------- Categorical Variable - Terrestrial Ecoregions --------
# read in raster biome
biome <- raster("C:/Data/niche-modelling/SDM_22/merremia/data/biome/wrld_biome.tif")
# crop the biome raster to the same extent as the soil- remove polar regions 
biome <- crop(biome, extent(slp), mask = T)
# resample biome layer to match extent, crs and resolution of elevation_wrld
biome <-  raster::projectRaster(biome, slp, method = "ngb")

# assign projection from slp layer to biome dataset
writeRaster(biome, filename="biome.tif", format ='GTiff', overwrite = TRUE)

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
install.packages('sdmpredictors')
library(sdmpredictors)

# run correlation co-efficient
p_mat <- pearson_correlation_matrix(stck, same_mask = FALSE)
View(p_mat) # remove one set of predictors that high correlation > 0.8 

# drop highly correlated predictors
stck = dropLayer(stck, c("bio1.tif", "bio6.tif", "bio17.tif", "world_aspect.tif", 
                         "preciptn.tif", "tempavg.tif", "tempmin.tif", "bio13.tif"))

# --------- GENERATE NATIVE RANGE PREDICTOR --------------------
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

# write the native range predictors to file
# change wd to native range data folder
setwd = "C:/Data/niche-modelling/SDM_22/merremia/30sec/nr_data"
writeRaster(stack(nr_predictor), names(nr_predictor), bylayer = TRUE,
            format = 'GTiff', overwrite = TRUE)
